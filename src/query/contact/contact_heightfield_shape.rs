use crate::bounding_volume::BoundingVolume;
use crate::math::{Pose, Real};
use crate::query::{Contact, QueryDispatcher};
use crate::shape::{HeightField, Shape};

/// Best contact between a heightfield and any other shape.
///
/// Tests the elements intersecting `shape2`'s Aabb (loosened by `prediction`) and returns
/// the deepest contact.
pub fn contact_heightfield_shape<D>(
    dispatcher: &D,
    pos12: &Pose,
    heightfield1: &HeightField,
    shape2: &dyn Shape,
    prediction: Real,
) -> Option<Contact>
where
    D: ?Sized + QueryDispatcher,
{
    let aabb2_1 = shape2.compute_aabb(pos12).loosened(prediction.max(0.0));
    let mut result = None::<Contact>;

    heightfield1.map_elements_in_local_aabb(&aabb2_1, &mut |_, elt1| {
        // The elements are already in the heightfield's local frame, so the sub-shape
        // contact needs no transform.
        if let Ok(Some(c)) = dispatcher.contact(pos12, elt1, shape2, prediction) {
            if result.is_none_or(|best| c.dist < best.dist) {
                result = Some(c);
            }
        }
    });

    result
}

/// Best contact between a shape and a heightfield.
pub fn contact_shape_heightfield<D>(
    dispatcher: &D,
    pos12: &Pose,
    shape1: &dyn Shape,
    heightfield2: &HeightField,
    prediction: Real,
) -> Option<Contact>
where
    D: ?Sized + QueryDispatcher,
{
    contact_heightfield_shape(
        dispatcher,
        &pos12.inverse(),
        heightfield2,
        shape1,
        prediction,
    )
    .map(|c| c.flipped())
}
