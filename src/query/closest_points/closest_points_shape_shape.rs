use crate::math::{Isometry, Real};
use crate::query::{ClosestPoints, DefaultQueryDispatcher, QueryDispatcher, Unsupported};
use crate::shape::Shape;

/// Computes the pair of closest points between two shapes.
///
/// Returns `None` if the objects are separated by a distance greater than `max_dist`.
pub fn closest_points(
    pos12: &Isometry<Real>,
    g1: &dyn Shape,
    g2: &dyn Shape,
    max_dist: Real,
) -> Result<ClosestPoints, Unsupported> {
    DefaultQueryDispatcher.closest_points(pos12, g1, g2, max_dist)
}
