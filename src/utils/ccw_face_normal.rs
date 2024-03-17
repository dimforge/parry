use crate::math::*;

/// Computes the direction pointing toward the right-hand-side of an oriented segment.
///
/// Returns `None` if the segment is degenerate.
#[inline]
#[cfg(feature = "dim2")]
pub fn ccw_face_normal(pts: [&Point; 2]) -> Option<UnitVector> {
    let ab = *pts[1] - *pts[0];
    let res = Vector::new(ab[1], -ab[0]);

    UnitVector::try_new(res, DEFAULT_EPSILON)
}

/// Computes the normal of a counter-clock-wise triangle.
///
/// Returns `None` if the triangle is degenerate.
#[inline]
#[cfg(feature = "dim3")]
pub fn ccw_face_normal(pts: [&Point; 3]) -> Option<UnitVector> {
    let ab = *pts[1] - *pts[0];
    let ac = *pts[2] - *pts[0];
    let res = ab.cross(ac);

    UnitVector::try_new(res, DEFAULT_EPSILON)
}
