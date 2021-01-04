use crate::math::{Isometry, Real};
use crate::query::{ClosestPoints, DefaultQueryDispatcher, QueryDispatcher, Unsupported};
use crate::shape::Shape;

/// Computes the pair of closest points between two shapes.
///
/// Returns `ClosestPoints::Disjoint` if the objects are separated by a distance greater than `max_dist`.
/// The result points in `ClosestPoints::WithinMargin` are expressed in world-space.
pub fn closest_points(
    pos1: &Isometry<Real>,
    g1: &dyn Shape,
    pos2: &Isometry<Real>,
    g2: &dyn Shape,
    max_dist: Real,
) -> Result<ClosestPoints, Unsupported> {
    let pos12 = pos1.inv_mul(pos2);
    DefaultQueryDispatcher
        .closest_points(&pos12, g1, g2, max_dist)
        .map(|res| res.transform_by(pos1, pos2))
}
