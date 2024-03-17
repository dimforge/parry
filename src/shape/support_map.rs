//! Traits for support mapping based shapes.

use crate::math::*;

/// Traits of convex shapes representable by a support mapping function.
///
/// # Parameters:
///   * V - type of the support mapping direction argument and of the returned point.
pub trait SupportMap {
    // Evaluates the support function of this shape.
    //
    // A support function is a function associating a vector to the shape point which maximizes
    // their dot product.
    fn local_support_point(&self, dir: &Vector) -> Point;

    /// Same as `self.local_support_point` except that `dir` is normalized.
    fn local_support_point_toward(&self, dir: &UnitVector) -> Point {
        self.local_support_point(&dir.into_inner())
    }

    // Evaluates the support function of this shape transformed by `transform`.
    //
    // A support function is a function associating a vector to the shape point which maximizes
    // their dot product.
    fn support_point(&self, transform: &Isometry, dir: &Vector) -> Point {
        let local_dir = transform.inverse_transform_vector(dir);
        transform.transform_point(&self.local_support_point(&local_dir))
    }

    /// Same as `self.support_point` except that `dir` is normalized.
    fn support_point_toward(&self, transform: &Isometry, dir: &UnitVector) -> Point {
        let local_dir = UnitVector::new_unchecked(transform.inverse_transform_vector(dir));
        transform.transform_point(&self.local_support_point_toward(&local_dir))
    }
}
