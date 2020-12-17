use na::Unit;

use crate::math::{Isometry, Real, Vector};
use crate::query::gjk::{self, DilatedShape, GJKResult, VoronoiSimplex};
use crate::query::{self, TOIStatus, TOI};
use crate::shape::SupportMap;
use num::Zero;

/// Time of impacts between two support-mapped shapes under translational movement.
pub fn time_of_impact_support_map_support_map<G1: ?Sized, G2: ?Sized>(
    pos12: &Isometry<Real>,
    vel12: &Vector<Real>,
    g1: &G1,
    g2: &G2,
    max_toi: Real,
    target_distance: Real,
) -> Option<TOI>
where
    G1: SupportMap,
    G2: SupportMap,
{
    if target_distance.is_zero() {
        gjk::directional_distance(pos12, g1, g2, &vel12, &mut VoronoiSimplex::new()).and_then(
            |(toi, normal1, witness1, witness2)| {
                if toi > max_toi {
                    None
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
    } else {
        let dilated1 = DilatedShape {
            shape: g1,
            radius: target_distance,
        };

        gjk::directional_distance(pos12, &dilated1, g2, &vel12, &mut VoronoiSimplex::new())
            .and_then(|(toi, normal1, witness1, witness2)| {
                if toi > max_toi {
                    None
                } else {
                    // This is mutable because, if the TOI is zero, we have to determine
                    // if the status isn't really a TOIStatus::Penetrating.
                    let mut status = TOIStatus::Converged;

                    if toi.is_zero() {
                        // The TOI is zero but we don't have valid witness points and normal
                        // yet because of we based our computations so far on the dilated shape.
                        // Therefore we need an extra step to retrieve the actual closest points, and
                        // determine if the actual shapes are penetrating or not.
                        // FIXME: all those computations are costly. Add a variant that returns only
                        // the TOI and does not computes the normal and witness points?
                        match query::details::closest_points_support_map_support_map_with_params(
                            pos12,
                            g1,
                            g2,
                            target_distance,
                            &mut VoronoiSimplex::new(),
                            Some(normal1),
                        ) {
                            GJKResult::ClosestPoints(pt1, pt2, _) => {
                                // Ok, we managed to compute the witness points.
                                let normal1 = Unit::new_normalize(pt2 - pt1);
                                return Some(TOI {
                                    toi,
                                    normal1,
                                    normal2: Unit::new_unchecked(
                                        pos12.inverse_transform_vector(&-normal1),
                                    ),
                                    witness1,
                                    witness2: pos12.inverse_transform_point(&pt2),
                                    status: TOIStatus::Converged,
                                });
                            }
                            GJKResult::NoIntersection(_) => {
                                // This should never happen.
                            }
                            GJKResult::Intersection => status = TOIStatus::Penetrating,
                            GJKResult::Proximity(_) => unreachable!(),
                        }
                    }

                    Some(TOI {
                        toi,
                        normal1: Unit::new_unchecked(normal1),
                        normal2: Unit::new_unchecked(pos12.inverse_transform_vector(&-normal1)),
                        witness1: witness1 - normal1 * target_distance,
                        witness2: pos12.inverse_transform_point(&witness2),
                        status,
                    })
                }
            })
    }
}
