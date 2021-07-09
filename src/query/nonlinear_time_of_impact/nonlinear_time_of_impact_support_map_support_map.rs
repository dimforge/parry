use na::{RealField, Unit};

use crate::math::{Point, Real, Vector};
use crate::query::{self, ClosestPoints, NonlinearRigidMotion, QueryDispatcher, TOIStatus, TOI};
use crate::shape::{Shape, SupportMap};
use crate::utils::WCross;

use crate::query::gjk::ConstantPoint;
use num::Bounded;

/// Enum specifying the behavior of TOI computation when there is a penetration at the starting time.
#[derive(Copy, Clone, Debug)]
pub enum NonlinearTOIMode {
    /// Stop TOI computation as soon as there is a penetration.
    StopAtPenetration,
    /// When there is a penetration, don't stop the TOI search if the relative velocity
    /// at the penetration points is negative (i.e. if the points are separating).
    DirectionalTOI {
        /// The sum of the `Shape::ccd_thickness` of both shapes involved in the TOI computation.
        sum_linear_thickness: Real,
        /// The max of the `Shape::ccd_angular_thickness` of both shapes involved in the TOI computation.
        max_angular_thickness: Real,
    },
}

impl NonlinearTOIMode {
    /// Initializes a directional TOI mode.
    ///
    /// With the "directional" TOI mode, the nonlinear TOI computation won't
    /// immediately stop if the shapes are already intersecting at `t = 0`.
    /// Instead, it will search for the first time where a contact between
    /// the shapes would result in a deeper penetration (with risk of tunnelling).
    /// This effectively checks the relative velocity of the shapes at their point
    /// of impact.
    pub fn directional_toi<S1, S2>(shape1: &S1, shape2: &S2) -> Self
    where
        S1: ?Sized + Shape,
        S2: ?Sized + Shape,
    {
        let sum_linear_thickness = shape1.ccd_thickness() + shape2.ccd_thickness();
        let max_angular_thickness = shape1
            .ccd_angular_thickness()
            .max(shape2.ccd_angular_thickness());

        NonlinearTOIMode::DirectionalTOI {
            sum_linear_thickness,
            max_angular_thickness,
        }
    }
}

/// Compute the time of first impact between two support-map shapes following
/// a nonlinear (with translations and rotations) motion.
pub fn nonlinear_time_of_impact_support_map_support_map<D, SM1, SM2>(
    dispatcher: &D,
    motion1: &NonlinearRigidMotion,
    sm1: &SM1,
    g1: &dyn Shape,
    motion2: &NonlinearRigidMotion,
    sm2: &SM2,
    g2: &dyn Shape,
    start_time: Real,
    end_time: Real,
    mode: NonlinearTOIMode,
) -> Option<TOI>
where
    D: ?Sized + QueryDispatcher,
    SM1: ?Sized + SupportMap,
    SM2: ?Sized + SupportMap,
{
    let sphere1 = g1.compute_local_bounding_sphere();
    let sphere2 = g2.compute_local_bounding_sphere();

    // Use the shape with the largest radius as the first shape.
    // This will give better convergence because everything will
    // be expressed in the local-space as that first shape (including
    // the separating axes used in the bisection).
    if sphere1.radius >= sphere2.radius {
        compute_toi(
            dispatcher, motion1, sm1, g1, motion2, sm2, g2, start_time, end_time, mode,
        )
    } else {
        compute_toi(
            dispatcher, motion2, sm2, g2, motion1, sm1, g1, start_time, end_time, mode,
        )
        .map(|toi| toi.swapped())
    }
}

/// Time of impacts between two support-mapped shapes under a rigid motion.
pub fn compute_toi<D, SM1, SM2>(
    dispatcher: &D,
    motion1: &NonlinearRigidMotion,
    sm1: &SM1,
    g1: &dyn Shape,
    motion2: &NonlinearRigidMotion,
    sm2: &SM2,
    g2: &dyn Shape,
    start_time: Real,
    end_time: Real,
    mode: NonlinearTOIMode,
) -> Option<TOI>
where
    D: ?Sized + QueryDispatcher,
    SM1: ?Sized + SupportMap,
    SM2: ?Sized + SupportMap,
{
    let mut prev_min_t = start_time;
    let abs_tol: Real = query::gjk::eps_tol();

    let mut result = TOI {
        toi: start_time,
        normal1: Vector::<Real>::x_axis(),
        normal2: Vector::<Real>::x_axis(),
        witness1: Point::<Real>::origin(),
        witness2: Point::<Real>::origin(),
        status: TOIStatus::Penetrating,
    };

    loop {
        let pos1 = motion1.position_at_time(result.toi);
        let pos2 = motion2.position_at_time(result.toi);
        let pos12 = pos1.inv_mul(&pos2);

        // TODO: use the _with_params version of the closest points query.
        match dispatcher
            .closest_points(&pos12, g1, g2, Real::max_value())
            .ok()?
        {
            ClosestPoints::Intersecting => {
                // println!(">> Intersecting.");
                if result.toi == start_time {
                    result.status = TOIStatus::Penetrating
                } else {
                    result.status = TOIStatus::Failed;
                }
                break;
            }
            ClosestPoints::WithinMargin(p1, p2) => {
                // println!(">> Within margin.");
                result.witness1 = p1;
                result.witness2 = p2;

                if let Some((normal1, dist)) =
                    Unit::try_new_and_get(pos12 * p2 - p1, crate::math::DEFAULT_EPSILON)
                {
                    // FIXME: do the "inverse transform unit vector" only when we are about to return.
                    result.normal1 = normal1;
                    result.normal2 = pos12.inverse_transform_unit_vector(&-normal1);

                    let curr_range = BisectionRange {
                        min_t: result.toi,
                        max_t: end_time,
                        curr_t: result.toi,
                    };

                    let (new_range, niter) =
                        bisect(dist, motion1, sm1, motion2, sm2, &normal1, curr_range);
                    // println!(
                    //     "Bisection result: {:?}, normal1: {:?}, normal2: {:?}",
                    //     new_range, result.normal1, result.normal2
                    // );
                    result.toi = new_range.curr_t;

                    if new_range.min_t - prev_min_t < abs_tol {
                        if new_range.max_t == end_time {
                            // Check the configuration at max_t to see if the object are not disjoint.
                            // NOTE: could we do this earlier, before the above loop?
                            // It feels like this could prevent catching some corner-cases like
                            // if one object is rotated by almost 180 degrees while the other is immobile.
                            let pos1 = motion1.position_at_time(new_range.max_t);
                            let pos2 = motion2.position_at_time(new_range.max_t);
                            let pos12 = pos1.inv_mul(&pos2);

                            let pt1 = sm1.local_support_point_toward(&normal1);
                            let pt2 = sm2.support_point_toward(&pos12, &-normal1);

                            if (pt2 - pt1).dot(&normal1) > 0.0 {
                                // We found an axis that separate both objects at the end configuration.
                                return None;
                            }
                        }

                        result.status = TOIStatus::Converged;
                        break;
                    }

                    prev_min_t = new_range.min_t;

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

    // In we started with a penetration, we need to compute a full contact manifold and
    // see if any of these contact points may result in tunnelling. If the is one, return
    // that TOI instead. That way, object moving tangentially on a surface (always keeping
    // a contact with it) won't report an useless impact.
    //
    // Note that this must be done here instead of outside of the `nonlinear_time_of_impact`
    // function so that this works properly with composite shapes.
    match mode {
        NonlinearTOIMode::DirectionalTOI {
            sum_linear_thickness,
            max_angular_thickness,
        } => {
            if (result.toi - start_time).abs() < 1.0e-5 {
                handle_penetration_at_start_time(
                    dispatcher,
                    motion1,
                    sm1,
                    g1,
                    motion2,
                    sm2,
                    g2,
                    start_time,
                    end_time,
                    sum_linear_thickness,
                    max_angular_thickness,
                )
            } else {
                Some(result)
            }
        }
        NonlinearTOIMode::StopAtPenetration => Some(result),
    }
}

fn handle_penetration_at_start_time<D, SM1, SM2>(
    dispatcher: &D,
    motion1: &NonlinearRigidMotion,
    sm1: &SM1,
    g1: &dyn Shape,
    motion2: &NonlinearRigidMotion,
    sm2: &SM2,
    g2: &dyn Shape,
    start_time: Real,
    end_time: Real,
    sum_linear_thickness: Real,
    max_angular_thickness: Real,
) -> Option<TOI>
where
    D: ?Sized + QueryDispatcher,
    SM1: ?Sized + SupportMap,
    SM2: ?Sized + SupportMap,
{
    // Because we are doing non-linear CCD, we need an iterative methode here.
    // First we need to check if the `toi = start_time` is legitimate, i.e.,
    // if tunnelling will happen if we don't clamp the motion.
    //
    // If the contact isn't "legitimate" (i.e. if we have a separating velocity),
    // then we need to check by how much the bodies can move without tunnelling.
    // With linear CCD it's easy: any linear movement is OK because we have
    // a separating velocity.
    //
    // With non-linear CCD it's more complicated because we need to take the
    // angular velocity into account, which will result in new contact points
    // that could tunnel (imagine for example 2D cuboid touching the ground at one
    // point. By rotating, its second end point will end up touching the ground too,
    // and we need to detect that so we don't permit a rotation larger than what's
    // needed for this second contact to happen).
    //
    // The iterative method here will iteratively check multiple rotation angles to
    // find new future contact points after some rotation; and check the relative
    // velocity at these future contact points.
    #[cfg(feature = "dim2")]
    let dangvel = (motion2.angvel - motion1.angvel).abs();
    #[cfg(feature = "dim3")]
    let dangvel = (motion2.angvel - motion1.angvel).norm();
    let inv_dangvel = crate::utils::inv(dangvel);
    let linear_increment = sum_linear_thickness;
    let angular_increment = Real::pi() - max_angular_thickness;

    let linear_time_increment =
        linear_increment * crate::utils::inv((motion2.linvel - motion1.linvel).norm());
    let angular_time_increment = angular_increment * inv_dangvel;
    let mut time_increment = angular_time_increment
        .min(linear_time_increment)
        // This is needed to avoid some tunnelling. But this is
        // kind of "brute force" so we should find something better.
        .min((end_time - start_time) / 10.0);

    // println!(
    //     "Lin time incr: {}, ang time incr: {}",
    //     linear_time_increment, angular_time_increment
    // );

    if time_increment == 0.0 {
        time_increment = end_time;
    }

    let mut next_time = start_time;

    // TODO: looping until we reach π sounds enough for most purposes.
    //       Is there a practical case where we need to loop until we reach 2π ?
    while next_time < end_time {
        // dbg!("A");
        let pos1_at_next_time = motion1.position_at_time(next_time);
        let pos2_at_next_time = motion2.position_at_time(next_time);
        let pos12_at_next_time = pos1_at_next_time.inv_mul(&pos2_at_next_time);
        let contact = dispatcher
            .contact(&pos12_at_next_time, g1, g2, Real::MAX)
            .ok()??;
        {
            // dbg!("C");

            // 1. Compute the relative velocity at that contact point.
            // 2. Check if this results in a potential tunnelling.
            // 3. Use bisection to adjust the TOI to the time where a pair
            //    of contact points potentially causing tunneling hit for the first time.
            let r1 = contact.point1 - motion1.local_center;
            let r2 = contact.point2 - motion2.local_center;
            let vel1 = motion1.linvel + motion1.angvel.gcross(pos1_at_next_time * r1);
            let vel2 = motion2.linvel + motion2.angvel.gcross(pos2_at_next_time * r2);
            let vel12 = vel2 - vel1;
            let normal_vel = -vel12.dot(&(pos1_at_next_time * contact.normal1));
            let ccd_threshold = if contact.dist <= 0.0 {
                sum_linear_thickness
            } else {
                contact.dist + sum_linear_thickness
            };

            // println!(
            //     "linvel: {:?}, angvel: {:?}, r2: {:?}, angpart: {:?}, vel2: {:?}",
            //     motion2.linvel,
            //     motion2.angvel,
            //     r2,
            //     motion2.angvel.gcross(pos2_at_next_time * r2),
            //     vel2,
            // );
            //
            // println!(
            //     "Found normal vel: {}, dist: {}, threshold: {}, if_value: {}, time: {}",
            //     normal_vel,
            //     contact.dist,
            //     ccd_threshold,
            //     normal_vel * (end_time - next_time),
            //     next_time
            // );

            if normal_vel * (end_time - next_time) > ccd_threshold {
                // dbg!("D1");

                let mut result = TOI {
                    toi: next_time,
                    witness1: contact.point1,
                    witness2: contact.point2,
                    normal1: contact.normal1,
                    normal2: contact.normal2,
                    status: TOIStatus::Converged,
                };

                if contact.dist > 0.0 {
                    // This is an acceptable impact. Now determine when
                    // the impacts happens exactly.
                    let curr_range = BisectionRange {
                        min_t: next_time,
                        max_t: end_time,
                        curr_t: next_time,
                    };
                    let (new_range, _) = bisect(
                        contact.dist,
                        motion1,
                        sm1,
                        motion2,
                        sm2,
                        &contact.normal1,
                        curr_range,
                    );

                    // TODO: the bisection isn't always enough here. We should check that we
                    //       still have a contact now. If not, we should run the loop from
                    //       nonlinear_time_of_impact_support_map_support_map again from this
                    //       point forward.
                    result.toi = new_range.curr_t;
                } else {
                    // dbg!("Bissecting points");
                    // This is an acceptable impact. Now determine when
                    // the impacts happens exactly.
                    let curr_range = BisectionRange {
                        min_t: start_time,
                        max_t: next_time,
                        curr_t: next_time,
                    };
                    let (new_range, _) = bisect(
                        contact.dist,
                        motion1,
                        &ConstantPoint(contact.point1),
                        motion2,
                        &ConstantPoint(contact.point2),
                        &contact.normal1,
                        curr_range,
                    );

                    // TODO: the bisection isn't always enough here. We should check that we
                    //       still have a contact now. If not, we should run the loop from
                    //       nonlinear_time_of_impact_support_map_support_map again from this
                    //       point forward.
                    result.toi = new_range.curr_t;
                }

                // println!("Fount new toi: {}", result.toi);
                return Some(result);
            }

            // dbg!("D2");
        }

        // If there is no angular velocity, we don't have to
        // continue because we can't rotate the object.
        if inv_dangvel == 0.0 {
            return None;
        }

        // dbg!("E");
        next_time += time_increment;
    }

    None
}

#[derive(Copy, Clone, Debug)]
struct BisectionRange {
    min_t: Real,
    curr_t: Real,
    max_t: Real,
}

fn bisect<SM1, SM2>(
    mut dist: Real,
    motion1: &NonlinearRigidMotion,
    sm1: &SM1,
    motion2: &NonlinearRigidMotion,
    sm2: &SM2,
    normal1: &Unit<Vector<Real>>,
    mut range: BisectionRange,
) -> (BisectionRange, usize)
where
    SM1: ?Sized + SupportMap,
    SM2: ?Sized + SupportMap,
{
    let abs_tol: Real = query::gjk::eps_tol();
    let rel_tol = abs_tol; // ComplexField::sqrt(abs_tol);
    let mut niter = 0;

    // Use the world-space normal so it doesn't move with the shapes.
    // This is necessary to reduce the risk of extracting a root that
    // is not the root happening at the smallest time.
    let pos1 = motion1.position_at_time(range.curr_t);
    let world_normal1 = pos1 * normal1;

    loop {
        // println!("Bisection dist: {}, range: {:?}", dist, range);
        // TODO: use the secant method too for finding the next iterate and converge more quickly.
        if dist < 0.0 {
            // Too close or penetration, go back in time.
            range.max_t = range.curr_t;
            range.curr_t = (range.min_t + range.curr_t) * 0.5;
        } else if dist > rel_tol {
            // Too far apart, go forward in time.
            range.min_t = range.curr_t;
            range.curr_t = (range.curr_t + range.max_t) * 0.5;
        } else {
            // Reached tolerance, break.
            // println!("Bisection, break on dist tolerance.");
            break;
        }

        if range.max_t - range.min_t < abs_tol {
            range.curr_t = range.max_t;
            // println!("Bisection, break on tiny range.");
            break;
        }

        let pos1 = motion1.position_at_time(range.curr_t);
        let pos2 = motion2.position_at_time(range.curr_t);
        let pos12 = pos1.inv_mul(&pos2);

        let normal1 = pos1.inverse_transform_unit_vector(&world_normal1);
        let pt1 = sm1.local_support_point_toward(&normal1);
        let pt2 = sm2.support_point_toward(&pos12, &-normal1);
        dist = pt2.coords.dot(&normal1) - pt1.coords.dot(&normal1);

        niter += 1;
    }

    (range, niter)
}
