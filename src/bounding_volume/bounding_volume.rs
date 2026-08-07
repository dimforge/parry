use crate::math::{Real, Vector};

/// Trait of bounding volumes.
///
/// Bounding volumes are coarse approximations of shapes. It usually have constant time
/// intersection, inclusion test. Two bounding volume must also be mergeable into a bigger bounding
/// volume.
pub trait BoundingVolume {
    // TODO: keep that ? What about non-spacial bounding volumes (e.g. bounding cones, curvature
    // bounds, etc.) ?
    /// Returns a point inside of this bounding volume. This is ideally its center.
    fn center(&self) -> Vector;

    /// Checks if this bounding volume intersect with another one.
    fn intersects(&self, _: &Self) -> bool;

    /// Checks if this bounding volume contains another one.
    fn contains(&self, _: &Self) -> bool;

    /// Merges this bounding volume with another one. The merge is done in-place.
    ///
    /// Due to floating-point rounding, the merged volume is not guaranteed to strictly
    /// [`contain`](Self::contains) both input volumes. Use [`loosened`](Self::loosened)
    /// with a small margin if strict containment is required.
    fn merge(&mut self, _: &Self);

    /// Merges this bounding volume with another one.
    ///
    /// Due to floating-point rounding, the merged volume is not guaranteed to strictly
    /// [`contain`](Self::contains) both input volumes. Use [`loosened`](Self::loosened)
    /// with a small margin if strict containment is required.
    fn merged(&self, _: &Self) -> Self;

    /// Enlarges this bounding volume.
    fn loosen(&mut self, _: Real);

    /// Creates a new, enlarged version, of this bounding volume.
    fn loosened(&self, _: Real) -> Self;

    /// Tighten this bounding volume.
    fn tighten(&mut self, _: Real);

    /// Creates a new, tightened version, of this bounding volume.
    fn tightened(&self, _: Real) -> Self;
}
