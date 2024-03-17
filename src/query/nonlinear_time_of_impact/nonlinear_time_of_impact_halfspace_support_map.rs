use na::RealField;

use crate::math::*;
use crate::query::{Ray, RayCast};
use crate::shape::HalfSpace;
use crate::shape::SupportMap;

/// Time Of Impact of a halfspace with a support-mapped shape under a rigid motion (translation + rotation).
pub fn nonlinear_time_of_impact_halfspace_support_map<G: ?Sized>(
    pos12: &Isometry,
    vel12: &Vector,
    halfspace: &HalfSpace,
    other: &G,
) -> Option
where
    G: SupportMap,
{
    /*
    let halfspace_normal = mhalfspace * halfspace.normal();
    let mut curr = 0.0;


    loop {
        let curr_mother = mother.advance(dvel);
        let closest_point = other.support_point(curr_mother, &-halfspace_normal);

    }

    let vel = *vel_other - *vel_halfspace;
    let closest_point = other.support_point(mother, &-halfspace_normal);

    halfspace.cast_ray(mhalfspace, &Ray::new(closest_point, vel), true)
    */
    unimplemented!()
}

/// Time Of Impact of a halfspace with a support-mapped shape under a rigid motion (translation + rotation).
pub fn nonlinear_time_of_impact_support_map_halfspace<G: ?Sized>(
    pos12: &Isometry,
    vel12: &Vector,
    other: &G,
    halfspace: &HalfSpace,
) -> Option
where
    G: SupportMap,
{
    nonlinear_time_of_impact_halfspace_support_map(pos12, vel12, halfspace, other)
}
