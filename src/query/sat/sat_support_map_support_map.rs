use crate::math::*;
use crate::shape::SupportMap;

/// Computes the separation along the given direction,
/// between two convex shapes implementing the `SupportMap` trait.
#[allow(dead_code)]
pub fn support_map_support_map_compute_separation(
    sm1: &impl SupportMap,
    sm2: &impl SupportMap,
    pos12: &Isometry,
    dir1: &UnitVector,
) -> Real {
    let p1 = sm1.local_support_point_toward(dir1);
    let p2 = sm2.support_point_toward(pos12, &-*dir1);
    (p2 - p1).dot(dir1.into_inner())
}
