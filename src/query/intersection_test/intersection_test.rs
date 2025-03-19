use crate::math::Isometry;
use crate::query::{DefaultQueryDispatcher, QueryDispatcher, Unsupported};
use crate::shape::Shape;

/// Tests whether two shapes are intersecting.
pub fn intersection_test(
    pos1: &Isometry,
    g1: &dyn Shape,
    pos2: &Isometry,
    g2: &dyn Shape,
) -> Result<bool, Unsupported> {
    let pos12 = pos1.inv_mul(pos2);
    DefaultQueryDispatcher.intersection_test(&pos12, g1, g2)
}
