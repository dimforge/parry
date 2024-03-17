use crate::math::*;
use crate::utils;
use na::{self, ComplexField};

/// Computes the bounding sphere of a set of point, given its center.
// FIXME: return a bounding sphere?
#[inline]
pub fn point_cloud_bounding_sphere_with_center(pts: &[Point], center: Point) -> (Point, Real) {
    let mut sqradius = 0.0;

    for pt in pts.iter() {
        let distance_squared = distance_squared(pt, center);

        if distance_squared > sqradius {
            sqradius = distance_squared
        }
    }

    (center, ComplexField::sqrt(sqradius))
}

/// Computes a bounding sphere of the specified set of point.
// FIXME: return a bounding sphere?
#[inline]
pub fn point_cloud_bounding_sphere(pts: &[Point]) -> (Point, Real) {
    point_cloud_bounding_sphere_with_center(pts, utils::center(pts))
}
