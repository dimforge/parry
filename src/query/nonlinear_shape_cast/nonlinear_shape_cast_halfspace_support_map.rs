use na::RealField;

use crate::math::{Isometry, Real, Vector};
use crate::query::{Ray, RayCast};
use crate::shape::HalfSpace;
use crate::shape::SupportMap;

/// Time Of Impact of a halfspace with a support-mapped shape under a rigid motion (translation + rotation).
pub fn cast_shapes_nonlinear_halfspace_support_map<G: ?Sized + SupportMap>(
    pos12: &Isometry<Real>,
    vel12: &Vector<Real>,
    halfspace: &HalfSpace,
    other: &G,
) -> Option {
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
pub fn cast_shapes_nonlinear_support_map_halfspace<G: ?Sized + SupportMap>(
    pos12: &Isometry<Real>,
    vel12: &Vector<Real>,
    other: &G,
    halfspace: &HalfSpace,
) -> Option {
    cast_shapes_nonlinear_halfspace_support_map(pos12, vel12, halfspace, other)
}
