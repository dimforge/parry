use crate::math::{Point, Real};
use crate::utils;
use na;

/// Computes the bounding sphere of a set of point, given its center.
// FIXME: return a bounding sphere?
#[inline]
pub fn point_cloud_bounding_sphere_with_center(
    pts: &[Point<Real>],
    center: Point<Real>,
) -> (Point<Real>, Real) {
    let mut sqradius = na::zero::<Real>();

    for pt in pts.iter() {
        let distance_squared = na::distance_squared(pt, &center);

        if distance_squared > sqradius {
            sqradius = distance_squared
        }
    }

    (center, sqradius.sqrt())
}

/// Computes a bounding sphere of the specified set of point.
// FIXME: return a bounding sphere?
#[inline]
pub fn point_cloud_bounding_sphere(pts: &[Point<Real>]) -> (Point<Real>, Real) {
    point_cloud_bounding_sphere_with_center(pts, utils::center(pts))
}
