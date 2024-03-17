use crate::math::*;

use crate::query::{DefaultQueryDispatcher, QueryDispatcher, Unsupported};
use crate::shape::Shape;

/// Computes the minimum distance separating two shapes.
///
/// Returns `0.0` if the objects are touching or penetrating.
pub fn distance(
    pos1: &Isometry,
    g1: &dyn Shape,
    pos2: &Isometry,
    g2: &dyn Shape,
) -> Result<Real, Unsupported> {
    let pos12 = pos1.inv_mul(pos2);
    DefaultQueryDispatcher.distance(&pos12, g1, g2)
}
