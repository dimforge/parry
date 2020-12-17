use crate::math::{Isometry, Real};
use crate::shape::HalfSpace;
use crate::shape::SupportMap;

/// Intersection test between a halfspace and a support-mapped shape (Cuboid, ConvexHull, etc.)
pub fn intersection_test_halfspace_support_map<G: ?Sized + SupportMap>(
    pos12: &Isometry<Real>,
    halfspace: &HalfSpace,
    other: &G,
) -> bool {
    let deepest = other.support_point_toward(pos12, &-halfspace.normal);
    halfspace.normal.dot(&deepest.coords) <= 0.0
}

/// Intersection test between a support-mapped shape (Cuboid, ConvexHull, etc.) and a halfspace.
pub fn intersection_test_support_map_halfspace<G: ?Sized + SupportMap>(
    pos12: &Isometry<Real>,
    other: &G,
    halfspace: &HalfSpace,
) -> bool {
    intersection_test_halfspace_support_map(&pos12.inverse(), halfspace, other)
}
