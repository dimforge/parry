use na::Unit;

use crate::math::{Isometry, Real, Vector};
use crate::query::details;
use crate::query::gjk::{self, VoronoiSimplex};
use crate::query::{TOIStatus, TOI};
use crate::shape::SupportMap;
use num::Zero;

/// Time of impacts between two support-mapped shapes under translational movement.
pub fn time_of_impact_support_map_support_map<G1: ?Sized, G2: ?Sized>(
    pos12: &Isometry<Real>,
    vel12: &Vector<Real>,
    g1: &G1,
    g2: &G2,
    max_toi: Real,
    stop_at_penetration: bool,
) -> Option<TOI>
where
    G1: SupportMap,
    G2: SupportMap,
{
    gjk::directional_distance(pos12, g1, g2, vel12, &mut VoronoiSimplex::new()).and_then(
        |(toi, normal1, witness1, witness2)| {
            if toi > max_toi {
                None
            } else if !stop_at_penetration && toi < 1.0e-5 {
                let contact = details::contact_support_map_support_map(pos12, g1, g2, Real::MAX)?;
                let normal_vel = contact.normal1.dot(vel12);

                if normal_vel >= 0.0 {
                    None
                } else {
                    Some(TOI {
                        toi,
                        normal1: contact.normal1,
                        normal2: contact.normal2,
                        witness1: contact.point1,
                        witness2: contact.point2,
                        status: TOIStatus::Penetrating,
                    })
                }
            } else {
                Some(TOI {
                    toi,
                    normal1: Unit::new_unchecked(normal1),
                    normal2: Unit::new_unchecked(pos12.inverse_transform_vector(&-normal1)),
                    witness1,
                    witness2: pos12.inverse_transform_point(&witness2),
                    status: if toi.is_zero() {
                        TOIStatus::Penetrating
                    } else {
                        TOIStatus::Converged
                    },
                })
            }
        },
    )
}
