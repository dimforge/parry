use crate::math::{Isometry, Real};
use crate::query::Contact;
use crate::shape::{HalfSpace, SupportMap};

/// Contact between a halfspace and a support-mapped shape (Cuboid, ConvexHull, etc.)
pub fn contact_halfspace_support_map<G: ?Sized + SupportMap>(
    pos12: &Isometry<Real>,
    halfspace: &HalfSpace,
    other: &G,
    prediction: Real,
) -> Option<Contact> {
    let deepest = other.support_point_toward(&pos12, &-halfspace.normal);
    let distance = halfspace.normal.dot(&deepest.coords);

    if distance <= prediction {
        let point1 = deepest - halfspace.normal.into_inner() * distance;
        let point2 = pos12.inverse_transform_point(&deepest);
        let normal2 = pos12.inverse_transform_unit_vector(&-halfspace.normal);

        Some(Contact::new(
            point1,
            point2,
            halfspace.normal,
            normal2,
            distance,
        ))
    } else {
        None
    }
}

/// Contact between a support-mapped shape (Cuboid, ConvexHull, etc.) and a halfspace.
pub fn contact_support_map_halfspace<G: ?Sized + SupportMap>(
    pos12: &Isometry<Real>,
    other: &G,
    halfspace: &HalfSpace,
    prediction: Real,
) -> Option<Contact> {
    contact_halfspace_support_map(pos12, halfspace, other, prediction).map(|c| c.flipped())
}
