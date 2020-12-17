use na::Unit;

use crate::math::{Isometry, Point, Real, Vector};
use crate::motion::RigidMotion;
use crate::query::{self, ClosestPoints, TOIStatus, TOI};
use crate::shape::SupportMap;

use num::Bounded;

/// Time of impacts between two support-mapped shapes under a rigid motion.
pub fn nonlinear_time_of_impact_support_map_support_map<G1: ?Sized, G2: ?Sized>(
    motion12: &(impl RigidMotion + ?Sized),
    g1: &G1,
    g2: &G2,
    max_toi: Real,
    target_distance: Real,
) -> Option<TOI>
where
    G1: SupportMap,
    G2: SupportMap,
{
    nonlinear_time_of_impact_support_map_support_map_with_closest_points_function(
        motion12,
        g1,
        g2,
        max_toi,
        target_distance,
        query::details::closest_points_support_map_support_map,
    )
}

/// Time of impacts between two support-mapped shapes under a rigid motion.
///
/// You probably want to use `query::details::nonlinear_time_of_impact_support_map_support_map` instead of this one.
/// The distance function between the two shapes must be given.
pub fn nonlinear_time_of_impact_support_map_support_map_with_closest_points_function<
    G1: ?Sized,
    G2: ?Sized,
>(
    motion12: &(impl RigidMotion + ?Sized),
    g1: &G1,
    g2: &G2,
    max_toi: Real,
    target_distance: Real,
    closest_points: impl Fn(&Isometry<Real>, &G1, &G2, Real) -> ClosestPoints,
) -> Option<TOI>
where
    G1: SupportMap,
    G2: SupportMap,
{
    let _0_5: Real = 0.5;
    let mut min_t = na::zero::<Real>();
    let mut prev_min_t = min_t;
    let abs_tol: Real = query::gjk::eps_tol();
    let rel_tol = abs_tol.sqrt();
    let mut result = TOI {
        toi: na::zero::<Real>(),
        normal1: Vector::<Real>::x_axis(),
        normal2: Vector::<Real>::x_axis(),
        witness1: Point::<Real>::origin(),
        witness2: Point::<Real>::origin(),
        status: TOIStatus::Penetrating,
    };

    loop {
        let pos12 = motion12.position_at_time(result.toi);

        // FIXME: use the _with_params version of the closest points query.
        match closest_points(&pos12, g1, g2, Real::max_value()) {
            ClosestPoints::Intersecting => {
                if result.toi == na::zero::<Real>() {
                    result.status = TOIStatus::Penetrating
                } else {
                    result.status = TOIStatus::Failed;
                }
                break;
            }
            ClosestPoints::WithinMargin(p1, p2) => {
                result.witness1 = p1;
                result.witness2 = p2;

                if let Some((dir, mut dist)) =
                    Unit::try_new_and_get(pos12 * p2 - p1, crate::math::DEFAULT_EPSILON)
                {
                    // FIXME: do the "inverse transform unit vector" only when we are about to return.
                    result.normal1 = dir;
                    result.normal2 = pos12.inverse_transform_unit_vector(&-dir);

                    let mut niter = 0;
                    min_t = result.toi;
                    let mut max_t = max_toi;
                    let min_target_distance = (target_distance - rel_tol).max(na::zero::<Real>());
                    let max_target_distance = target_distance + rel_tol;

                    loop {
                        // FIXME: use the secant method too for finding the next iterate.
                        if dist < min_target_distance {
                            // Too close or penetration, go back in time.
                            max_t = result.toi;
                            result.toi = (min_t + result.toi) * _0_5;
                        } else if dist > max_target_distance {
                            // Too far apart, go forward in time.
                            min_t = result.toi;
                            result.toi = (result.toi + max_t) * _0_5;
                        } else {
                            // Reached tolerance, break.
                            break;
                        }

                        if max_t - min_t < abs_tol {
                            result.toi = min_t;
                            break;
                        }

                        let pos12 = motion12.position_at_time(result.toi);
                        let pt1 = g1.local_support_point_toward(&dir);
                        let pt2 = g2.support_point_toward(&pos12, &-dir);

                        dist = (pt2 - pt1).dot(&dir);

                        niter += 1;
                    }

                    if min_t - prev_min_t < abs_tol {
                        if max_t == max_toi {
                            let pos12 = motion12.position_at_time(max_t);

                            // Check the configuration at max_t to see if the object are not disjoint.
                            // NOTE: could we do this earlier, before the above loop?
                            // It feels like this could prevent catching some corner-cases like
                            // if one object is rotated by almost 180 degrees while the other is immobile.
                            let pt1 = g1.local_support_point_toward(&dir);
                            let pt2 = g2.support_point_toward(&pos12, &-dir);
                            if (pt2 - pt1).dot(&dir) > target_distance {
                                // We found an axis that separate both objects at the end configuration.
                                return None;
                            }
                        }

                        result.status = TOIStatus::Converged;
                        break;
                    }

                    prev_min_t = min_t;

                    if niter == 0 {
                        result.status = TOIStatus::Converged;
                        break;
                    }
                } else {
                    result.status = TOIStatus::Failed;
                    break;
                }
            }
            ClosestPoints::Disjoint => unreachable!(),
        }
    }

    Some(result)
}
