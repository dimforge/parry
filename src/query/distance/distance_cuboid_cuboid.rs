use crate::math::{Isometry, Real, real};
use crate::query::ClosestPoints;
use crate::shape::Cuboid;

/// Distance between two cuboids.
#[inline]
pub fn distance_cuboid_cuboid(pos12: &Isometry<Real>, cuboid1: &Cuboid, cuboid2: &Cuboid) -> Real {
    match crate::query::details::closest_points_cuboid_cuboid(pos12, cuboid1, cuboid2, Real::MAX) {
        ClosestPoints::WithinMargin(p1, p2) => na::distance(&p1, &(pos12 * p2)),
        _ => real!(0.0),
    }
}
