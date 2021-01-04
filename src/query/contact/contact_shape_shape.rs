use crate::math::{Isometry, Real};
use crate::query::{Contact, DefaultQueryDispatcher, QueryDispatcher, Unsupported};
use crate::shape::Shape;

/// Computes one pair of contact points point between two shapes.
///
/// Returns `None` if the objects are separated by a distance greater than `prediction`.
/// The result is given in world-space.
pub fn contact(
    pos1: &Isometry<Real>,
    g1: &dyn Shape,
    pos2: &Isometry<Real>,
    g2: &dyn Shape,
    prediction: Real,
) -> Result<Option<Contact>, Unsupported> {
    let pos12 = pos1.inv_mul(pos2);
    let mut result = DefaultQueryDispatcher.contact(&pos12, g1, g2, prediction);

    if let Ok(Some(contact)) = &mut result {
        contact.transform_by_mut(pos1, pos2);
    }

    result
}
