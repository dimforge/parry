//! Traits for support mapping based shapes.

use crate::math::{Isometry, Point, Real, Vector};
use na::Unit;

/// Traits of convex shapes representable by a support mapping function.
///
/// # Parameters:
///   * V - type of the support mapping direction argument and of the returned point.
pub trait SupportMap {
    // Evaluates the support function of this shape.
    //
    // A support function is a function associating a vector to the shape point which maximizes
    // their dot product.
    fn local_support_point(&self, dir: &Vector<Real>) -> Point<Real>;

    /// Same as `self.local_support_point` except that `dir` is normalized.
    fn local_support_point_toward(&self, dir: &Unit<Vector<Real>>) -> Point<Real> {
        self.local_support_point(dir.as_ref())
    }

    // Evaluates the support function of this shape transformed by `transform`.
    //
    // A support function is a function associating a vector to the shape point which maximizes
    // their dot product.
    fn support_point(&self, transform: &Isometry<Real>, dir: &Vector<Real>) -> Point<Real> {
        let local_dir = transform.inverse_transform_vector(dir);
        transform * self.local_support_point(&local_dir)
    }

    /// Same as `self.support_point` except that `dir` is normalized.
    fn support_point_toward(
        &self,
        transform: &Isometry<Real>,
        dir: &Unit<Vector<Real>>,
    ) -> Point<Real> {
        let local_dir = Unit::new_unchecked(transform.inverse_transform_vector(dir));
        transform * self.local_support_point_toward(&local_dir)
    }
}
