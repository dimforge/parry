use crate::math::{Point, Real};
use na;

/// Computes the geometric center (centroid) of a set of points.
///
/// The center is calculated by averaging all the point coordinates. This is also known as
/// the centroid or barycenter of the point cloud. All points are weighted equally.
///
/// # Arguments
///
/// * `pts` - A slice of points. Must contain at least one point.
///
/// # Returns
///
/// The geometric center as a `Point<Real>`.
///
/// # Panics
///
/// Panics if the input slice is empty.
///
/// # Examples
///
/// ## 2D Example
///
/// ```
/// # #[cfg(all(feature = "dim2", feature = "f32"))] {
/// use parry2d::utils::center;
/// use parry2d::math::Point;
///
/// let points = vec![
///     Point::new(0.0, 0.0),
///     Point::new(2.0, 0.0),
///     Point::new(2.0, 2.0),
///     Point::new(0.0, 2.0),
/// ];
///
/// let c = center(&points);
///
/// // The center of a square is at its middle
/// assert!((c.x - 1.0).abs() < 1e-6);
/// assert!((c.y - 1.0).abs() < 1e-6);
/// # }
/// ```
///
/// ## 3D Example
///
/// ```
/// # #[cfg(all(feature = "dim3", feature = "f32"))] {
/// use parry3d::utils::center;
/// use parry3d::math::Point;
///
/// let points = vec![
///     Point::new(0.0, 0.0, 0.0),
///     Point::new(4.0, 0.0, 0.0),
///     Point::new(0.0, 4.0, 0.0),
/// ];
///
/// let c = center(&points);
///
/// // The center of these three points
/// assert!((c.x - 4.0 / 3.0).abs() < 1e-6);
/// assert!((c.y - 4.0 / 3.0).abs() < 1e-6);
/// assert!(c.z.abs() < 1e-6);
/// # }
/// ```
///
/// ## Single Point
///
/// ```
/// # #[cfg(all(feature = "dim2", feature = "f32"))] {
/// use parry2d::utils::center;
/// use parry2d::math::Point;
///
/// let points = vec![Point::new(5.0, 10.0)];
/// let c = center(&points);
///
/// // The center of a single point is the point itself
/// assert_eq!(c, points[0]);
/// # }
/// ```
#[inline]
pub fn center(pts: &[Point<Real>]) -> Point<Real> {
    assert!(
        !pts.is_empty(),
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
