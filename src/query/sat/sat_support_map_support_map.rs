use crate::math::{Isometry, Real, Vector};
use crate::shape::SupportMap;
use na::Unit;

/// Computes the separation along the given direction,
/// between two convex shapes implementing the `SupportMap` trait.
#[allow(dead_code)]
pub fn support_map_support_map_compute_separation(
    sm1: &impl SupportMap,
    sm2: &impl SupportMap,
    pos12: &Isometry<Real>,
    dir1: &Unit<Vector<Real>>,
) -> Real {
    let p1 = sm1.local_support_point_toward(dir1);
    let p2 = sm2.support_point_toward(pos12, &-*dir1);
    (p2 - p1).dot(dir1)
}
