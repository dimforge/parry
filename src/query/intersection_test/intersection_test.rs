use crate::math::{Isometry, Real};
use crate::query::{DefaultQueryDispatcher, QueryDispatcher, Unsupported};
use crate::shape::Shape;

/// Tests whether two shapes are intersecting.
pub fn intersection_test(
    pos12: &Isometry<Real>,
    g1: &dyn Shape,
    g2: &dyn Shape,
) -> Result<bool, Unsupported> {
    DefaultQueryDispatcher.intersection_test(pos12, g1, g2)
}
