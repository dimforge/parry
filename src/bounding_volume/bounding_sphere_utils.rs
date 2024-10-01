use crate::math::{Point, Real};
use crate::utils;
use na::{self, ComplexField};

use super::BoundingSphere;

/// Computes the bounding sphere of a set of point, given its center.
#[inline]
pub fn point_cloud_bounding_sphere_with_center(
    pts: &[Point<Real>],
    center: Point<Real>,
) -> BoundingSphere {
    let mut sqradius = 0.0;

    for pt in pts.iter() {
        let distance_squared = na::distance_squared(pt, &center);

        if distance_squared > sqradius {
            sqradius = distance_squared
        }
    }
    BoundingSphere::new(center, ComplexField::sqrt(sqradius))
}

/// Computes a bounding sphere of the specified set of point.
#[inline]
pub fn point_cloud_bounding_sphere(pts: &[Point<Real>]) -> BoundingSphere {
    point_cloud_bounding_sphere_with_center(pts, utils::center(pts))
}
