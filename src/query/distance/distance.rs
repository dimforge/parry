use crate::math::{Isometry, Real};

use crate::query::{DefaultQueryDispatcher, QueryDispatcher, Unsupported};
use crate::shape::Shape;

/// Computes the minimum distance separating two shapes.
///
/// Returns `0.0` if the objects are touching or penetrating.
pub fn distance(
    pos12: &Isometry<Real>,
    g1: &dyn Shape,
    g2: &dyn Shape,
) -> Result<Real, Unsupported> {
    DefaultQueryDispatcher.distance(pos12, g1, g2)
}
