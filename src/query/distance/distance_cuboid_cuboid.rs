use crate::math::*;
use crate::query::ClosestPoints;
use crate::shape::Cuboid;

/// Distance between two cuboids.
#[inline]
pub fn distance_cuboid_cuboid(pos12: &Isometry, cuboid1: &Cuboid, cuboid2: &Cuboid) -> Real {
    match crate::query::details::closest_points_cuboid_cuboid(pos12, cuboid1, cuboid2, Real::MAX) {
        ClosestPoints::WithinMargin(p1, p2) => distance(p1, pos12.transform_point(&p2)),
        _ => 0.0,
    }
}
