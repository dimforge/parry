use crate::math::{Isometry, Real};
use crate::query::{Contact, DefaultQueryDispatcher, QueryDispatcher, Unsupported};
use crate::shape::Shape;

/// Computes one pair of contact points point between two shapes.
///
/// Returns `None` if the objects are separated by a distance greater than `prediction`.
pub fn contact(
    pos12: &Isometry<Real>,
    g1: &dyn Shape,
    g2: &dyn Shape,
    prediction: Real,
) -> Result<Option<Contact>, Unsupported> {
    DefaultQueryDispatcher.contact(pos12, g1, g2, prediction)
}
