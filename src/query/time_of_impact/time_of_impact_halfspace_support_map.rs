use crate::math::{Isometry, Real, Vector};
use crate::query::{Ray, RayCast, TOIStatus, TOI};
use crate::shape::{HalfSpace, SupportMap};

/// Time Of Impact of a halfspace with a support-mapped shape under translational movement.
pub fn time_of_impact_halfspace_support_map<G: ?Sized>(
    pos12: &Isometry<Real>,
    vel12: &Vector<Real>,
    halfspace: &HalfSpace,
    other: &G,
    max_toi: Real,
) -> Option<TOI>
where
    G: SupportMap,
{
    // FIXME: add method to get only the local support point.
    // This would avoid the `inverse_transform_point` later.
    let support_point = other.support_point(pos12, &-halfspace.normal);
    let closest_point = support_point;
    let ray = Ray::new(closest_point, *vel12);

    if let Some(toi) = halfspace.cast_local_ray(&ray, max_toi, true) {
        if toi > max_toi {
            return None;
        }

        let status;
        let witness2 = support_point;
        let mut witness1 = ray.point_at(toi);

        if support_point.coords.dot(&halfspace.normal) < 0.0 {
            status = TOIStatus::Penetrating
        } else {
            // Project the witness point to the halfspace.
            // Note that witness2 is already in the halfspace's local-space.
            witness1 = witness1 - *halfspace.normal * witness1.coords.dot(&halfspace.normal);
            status = TOIStatus::Converged
        }

        Some(TOI {
            toi,
            normal1: halfspace.normal,
            normal2: pos12.inverse_transform_unit_vector(&-halfspace.normal),
            witness1,
            witness2: pos12.inverse_transform_point(&witness2),
            status,
        })
    } else {
        None
    }
}

/// Time Of Impact of a halfspace with a support-mapped shape under translational movement.
pub fn time_of_impact_support_map_halfspace<G: ?Sized>(
    pos12: &Isometry<Real>,
    vel12: &Vector<Real>,
    other: &G,
    halfspace: &HalfSpace,
    max_toi: Real,
) -> Option<TOI>
where
    G: SupportMap,
{
    time_of_impact_halfspace_support_map(
        &pos12.inverse(),
        &-pos12.inverse_transform_vector(&vel12),
        halfspace,
        other,
        max_toi,
    )
    .map(|toi| toi.swapped())
}
