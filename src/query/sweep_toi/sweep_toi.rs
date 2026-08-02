//! Conservative-advancement time of impact on endpoint-interpolated sweeps.
//!
//! A GJK distance query (with a warm-started simplex cache) provides a separating axis, and a
//! push-back loop with a hybrid bisection/false-position root finder advances the earliest
//! time at which the shapes reach the target separation.

use super::proxy_distance::{proxy_distance, SimplexCache};
use super::separation::SeparationFunction;
use super::sweep::{rotate_vec, Sweep};
use super::toi_proxy::ToiProxy;
use crate::math::{Real, Vector};

/// The outcome of a [`sweep_time_of_impact`] computation.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum SweepToiStatus {
    /// The shapes were overlapped at the start time; continuous collision gave up
    /// (`fraction == 0`).
    Overlapped,
    /// The shapes reach the target separation at `fraction`.
    Hit,
    /// The shapes never come within the target separation over the sweep
    /// (`fraction == max_fraction`).
    Separated,
    /// The root finder failed to converge; `fraction` holds the last safe time.
    Failed,
}

/// The result of a [`sweep_time_of_impact`] computation.
#[derive(Copy, Clone, Debug)]
pub struct SweepToiOutput {
    /// How the computation terminated.
    pub status: SweepToiStatus,
    /// The normalized time of impact in `[0, max_fraction]`.
    pub fraction: Real,
    /// The averaged world-space hit point (only meaningful on `Hit`/`Failed`).
    pub point: Vector,
    /// The world-space separating axis at the hit time, pointing from the first shape to the
    /// second (only meaningful on `Hit`/`Failed`).
    pub normal: Vector,
}

/// Computes the time of impact between two proxies following endpoint-interpolated sweeps.
///
/// The shapes are advanced to the earliest time in `[0, max_fraction]` at which their
/// core-shape separation drops to `max(linear_slop, radius_a + radius_b - linear_slop)`.
/// `linear_slop` is typically 0.005 length units.
pub fn sweep_time_of_impact(
    proxy_a: &ToiProxy,
    sweep_a: &Sweep,
    proxy_b: &ToiProxy,
    sweep_b: &Sweep,
    max_fraction: Real,
    linear_slop: Real,
) -> SweepToiOutput {
    let mut output = SweepToiOutput {
        status: SweepToiStatus::Separated,
        fraction: max_fraction,
        point: Vector::ZERO,
        normal: Vector::ZERO,
    };

    // Shift to the first sweep’s start center for better floating-point accuracy.
    let origin = sweep_a.c1;
    let sweep_a = sweep_a.shifted(origin);
    let sweep_b = sweep_b.shifted(origin);

    #[cfg(feature = "dim2")]
    let max_push_back_iterations = 8; // Maximum polygon vertex count.
    #[cfg(feature = "dim3")]
    let max_push_back_iterations = proxy_a.points().len() + proxy_b.points().len();

    #[cfg(feature = "dim2")]
    const MAX_DISTANCE_ITERATIONS: u32 = 20;
    #[cfg(feature = "dim3")]
    const MAX_DISTANCE_ITERATIONS: u32 = 25;

    let t_max = max_fraction;

    // Set up target distance and tolerance.
    let total_radius = proxy_a.radius + proxy_b.radius;
    let target = linear_slop.max(total_radius - linear_slop);
    let tolerance = 0.25 * linear_slop;
    debug_assert!(target > tolerance);

    let mut t1 = 0.0;
    let mut distance_iterations = 0u32;
    let mut cache = SimplexCache::default();

    // The outer loop progressively attempts to compute new separating axes.
    // This loop terminates when an axis is repeated (no progress is made).
    loop {
        // Get the distance between shapes. We can also use the results to get a separating
        // axis.
        let xf_a = sweep_a.transform_at(t1);
        let xf_b = sweep_b.transform_at(t1);
        let pos12 = xf_a.inv_mul(&xf_b);
        let distance_output = proxy_distance(&pos12, proxy_a, proxy_b, false, &mut cache);

        // The distance query runs in frame A, project the witness data back to the (shifted)
        // world.
        let world_normal = rotate_vec(&xf_a.rotation, distance_output.normal);
        let world_point_a = xf_a.transform_point(distance_output.point_a);
        let world_point_b = xf_a.transform_point(distance_output.point_b);

        distance_iterations += 1;

        let averaged_hit_point = || {
            let pa = world_point_a + proxy_a.radius * world_normal;
            let pb = world_point_b - proxy_b.radius * world_normal;
            0.5 * (pa + pb) + origin
        };

        // If the shapes are overlapped, we give up on continuous collision.
        if distance_output.distance <= 0.0 {
            output.status = SweepToiStatus::Overlapped;
            output.fraction = 0.0;
            break;
        }

        if distance_output.distance <= target + tolerance {
            // Success!
            output.status = SweepToiStatus::Hit;
            output.point = averaged_hit_point();
            output.normal = world_normal;
            output.fraction = t1;
            break;
        }

        // In 3D, check for slow progress before running the push-back loop…
        #[cfg(feature = "dim3")]
        if distance_iterations == MAX_DISTANCE_ITERATIONS {
            // Progress too slow. This can happen when a capsule rotates around a triangle
            // vertex.
            output.status = SweepToiStatus::Failed;
            output.fraction = t1;
            output.point = averaged_hit_point();
            output.normal = world_normal;
            break;
        }

        // Initialize the separating axis.
        let mut fcn = SeparationFunction::new(
            &cache,
            proxy_a,
            &sweep_a,
            proxy_b,
            &sweep_b,
            world_normal,
            t1,
        );

        // Compute the TOI on the separating axis. We do this by successively resolving the
        // deepest point. This loop is bounded by the number of vertices.
        let mut done = false;
        let mut t2 = t_max;
        let mut push_back_iterations = 0;
        loop {
            // Find the deepest point at t2. Store the witness point indices.
            let (mut s2, index_a, index_b) = fcn.find_min_separation(t2);

            // Is the final configuration separated?
            if s2 - target > tolerance {
                // Victory!
                output.status = SweepToiStatus::Separated;
                output.fraction = t_max;
                done = true;
                break;
            }

            // Has the separation reached tolerance?
            if s2 >= target - tolerance {
                // Advance the sweeps.
                t1 = t2;
                break;
            }

            // Compute the initial separation of the witness points.
            let mut s1 = fcn.evaluate(index_a, index_b, t1);

            // Check for initial overlap. This might happen if the root finder runs out of
            // iterations.
            if s1 < target - tolerance {
                output.status = SweepToiStatus::Failed;
                output.fraction = t1;
                done = true;
                break;
            }

            // Check for touching.
            if s1 <= target + tolerance {
                // Success! t1 should hold the TOI (could be 0.0).
                output.status = SweepToiStatus::Hit;
                output.point = averaged_hit_point();
                output.normal = world_normal;
                output.fraction = t1;
                done = true;
                break;
            }

            // Compute 1D root of: f(t) - target = 0.
            let mut root_iteration_count = 0;
            const MAX_ROOT_ITERATIONS: u32 = 50;
            let mut a1 = t1;
            let mut a2 = t2;
            loop {
                // Use a mix of false position and bisection.
                let t = if root_iteration_count & 1 == 1 {
                    // False position to improve convergence.
                    a1 + (target - s1) * (a2 - a1) / (s2 - s1)
                } else {
                    // Bisection to guarantee progress.
                    0.5 * (a1 + a2)
                };

                root_iteration_count += 1;

                let s = fcn.evaluate(index_a, index_b, t);

                // Has the separation reached tolerance?
                if (s - target).abs() <= tolerance {
                    // t2 holds a tentative value for t1.
                    t2 = t;
                    break;
                }

                // Ensure we continue to bracket the root.
                if s > target {
                    a1 = t;
                    s1 = s;
                } else {
                    a2 = t;
                    s2 = s;
                }

                if root_iteration_count == MAX_ROOT_ITERATIONS {
                    break;
                }
            }

            // Restart the inner loop if we have a failing edge case (3D edge/edge axis only).
            if root_iteration_count == MAX_ROOT_ITERATIONS - 1 && fcn.uses_edge_axis() {
                t2 = t_max;
                fcn.force_fixed_axis(t1);
            }

            push_back_iterations += 1;
            if push_back_iterations == max_push_back_iterations {
                break;
            }
        }

        if done {
            break;
        }

        // …while in 2D it is checked after the push-back loop, letting the last distance
        // iteration still resolve an impact.
        #[cfg(feature = "dim2")]
        if distance_iterations == MAX_DISTANCE_ITERATIONS {
            // Root finder got stuck. Semi-victory.
            output.status = SweepToiStatus::Failed;
            output.point = averaged_hit_point();
            output.normal = world_normal;
            output.fraction = t1;
            break;
        }
    }

    output
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::{Pose, Rotation};

    fn slop() -> Real {
        0.005
    }

    #[cfg(feature = "dim2")]
    fn pose(x: Real, y: Real) -> Pose {
        Pose::from_parts(Vector::new(x, y), Rotation::identity())
    }
    #[cfg(feature = "dim3")]
    fn pose(x: Real, y: Real) -> Pose {
        Pose::from_parts(Vector::new(x, y, 0.0), Rotation::IDENTITY)
    }

    fn ball_proxy(radius: Real) -> ToiProxy<'static> {
        ToiProxy::point(Vector::ZERO, radius)
    }

    fn cuboid_proxy(half_extent: Real) -> ToiProxy<'static> {
        #[cfg(feature = "dim2")]
        {
            ToiProxy::from_array(
                [
                    Vector::new(-half_extent, -half_extent),
                    Vector::new(half_extent, -half_extent),
                    Vector::new(half_extent, half_extent),
                    Vector::new(-half_extent, half_extent),
                ],
                0.0,
            )
        }
        #[cfg(feature = "dim3")]
        {
            let h = half_extent;
            ToiProxy::from_array(
                [
                    Vector::new(-h, -h, -h),
                    Vector::new(h, -h, -h),
                    Vector::new(h, h, -h),
                    Vector::new(-h, h, -h),
                    Vector::new(-h, -h, h),
                    Vector::new(h, -h, h),
                    Vector::new(h, h, h),
                    Vector::new(-h, h, h),
                ],
                0.0,
            )
        }
    }

    #[test]
    fn ball_hits_static_box() {
        // Static unit box at the origin; ball of radius 0.1 sweeping from x = -3 to x = +3.
        let wall = cuboid_proxy(0.5);
        let wall_sweep = Sweep::constant(&pose(0.0, 0.0), Vector::ZERO);
        let ball = ball_proxy(0.1);
        let ball_sweep = Sweep::from_poses(&pose(-3.0, 0.0), &pose(3.0, 0.0), Vector::ZERO);

        let result = sweep_time_of_impact(&wall, &wall_sweep, &ball, &ball_sweep, 1.0, slop());
        assert_eq!(result.status, SweepToiStatus::Hit);

        // The ball is stopped when its core (center) is `target = radius - slop` away from
        // the box surface: center_x = -0.5 - (0.1 - 0.005), fraction = (start - x) / travel.
        let expected = (3.0 - 0.5 - (0.1 - slop())) / 6.0;
        assert!(
            (result.fraction - expected).abs() < 0.001,
            "fraction {} vs expected {expected}",
            result.fraction
        );
    }

    #[test]
    fn ball_misses_box() {
        // Ball sweeping parallel to the box, far away.
        let wall = cuboid_proxy(0.5);
        let wall_sweep = Sweep::constant(&pose(0.0, 0.0), Vector::ZERO);
        let ball = ball_proxy(0.1);
        let ball_sweep = Sweep::from_poses(&pose(-3.0, 5.0), &pose(3.0, 5.0), Vector::ZERO);

        let result = sweep_time_of_impact(&wall, &wall_sweep, &ball, &ball_sweep, 1.0, slop());
        assert_eq!(result.status, SweepToiStatus::Separated);
        assert_eq!(result.fraction, 1.0);
    }

    #[test]
    fn overlapped_at_start_returns_fraction_zero() {
        let wall = cuboid_proxy(0.5);
        let wall_sweep = Sweep::constant(&pose(0.0, 0.0), Vector::ZERO);
        let ball = ball_proxy(0.1);
        let ball_sweep = Sweep::from_poses(&pose(0.0, 0.0), &pose(3.0, 0.0), Vector::ZERO);

        let result = sweep_time_of_impact(&wall, &wall_sweep, &ball, &ball_sweep, 1.0, slop());
        assert_eq!(result.status, SweepToiStatus::Overlapped);
        assert_eq!(result.fraction, 0.0);
    }

    #[test]
    fn box_vs_box_face_impact() {
        // Two unit boxes; one sweeps into the other along x.
        let a = cuboid_proxy(0.5);
        let a_sweep = Sweep::constant(&pose(0.0, 0.0), Vector::ZERO);
        let b = cuboid_proxy(0.5);
        let b_sweep = Sweep::from_poses(&pose(-4.0, 0.0), &pose(4.0, 0.0), Vector::ZERO);

        let result = sweep_time_of_impact(&a, &a_sweep, &b, &b_sweep, 1.0, slop());
        assert_eq!(result.status, SweepToiStatus::Hit);

        // Faces meet when the gap reaches target = slop: center_x = -(0.5 + 0.5 + slop).
        let expected = (4.0 - 1.0 - slop()) / 8.0;
        assert!(
            (result.fraction - expected).abs() < 0.001,
            "fraction {} vs expected {expected}",
            result.fraction
        );
    }

    #[test]
    fn rotating_bar_hits_ball() {
        // A long thin bar rotating 90° about the origin must catch a ball placed along its
        // swept arc even though the bar’s endpoint poses don’t overlap the ball.
        let half_length = 2.0;
        #[cfg(feature = "dim2")]
        let (bar, start, end) = (
            ToiProxy::from_array(
                [Vector::new(-half_length, 0.0), Vector::new(half_length, 0.0)],
                0.05,
            ),
            Pose::from_parts(Vector::ZERO, Rotation::identity()),
            Pose::from_parts(Vector::ZERO, Rotation::from_angle(core::f32::consts::FRAC_PI_2 as Real)),
        );
        #[cfg(feature = "dim3")]
        let (bar, start, end) = (
            ToiProxy::from_array(
                [
                    Vector::new(-half_length, 0.0, 0.0),
                    Vector::new(half_length, 0.0, 0.0),
                ],
                0.05,
            ),
            Pose::from_parts(Vector::ZERO, Rotation::IDENTITY),
            Pose::from_parts(
                Vector::ZERO,
                Rotation::from_rotation_z(core::f32::consts::FRAC_PI_2 as Real),
            ),
        );

        let bar_sweep = Sweep::from_poses(&start, &end, Vector::ZERO);
        let ball = ball_proxy(0.1);
        // Place the ball at 45° on the arc of radius 1.5.
        let d = 1.5 * (0.5_f64.sqrt() as Real);
        let ball_sweep = Sweep::constant(&pose(d, d), Vector::ZERO);

        let result = sweep_time_of_impact(&bar, &bar_sweep, &ball, &ball_sweep, 1.0, slop());
        assert_eq!(result.status, SweepToiStatus::Hit);
        // The bar reaches 45° at t = 0.5; it should hit slightly before.
        assert!(
            result.fraction > 0.3 && result.fraction < 0.5,
            "fraction {}",
            result.fraction
        );
    }
}
