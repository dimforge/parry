use crate::math::{Point, Real};
use na;

/// Computes the center of a set of point.
#[inline]
pub fn center(pts: &[Point<Real>]) -> Point<Real> {
    assert!(
        pts.len() >= 1,
        "Cannot compute the center of less than 1 point."
    );

    let denom: Real = na::convert::<f64, Real>(1.0 / (pts.len() as f64));

    let mut piter = pts.iter();
    let mut res = *piter.next().unwrap() * denom;

    for pt in piter {
        res += pt.coords * denom;
    }

    res
}
