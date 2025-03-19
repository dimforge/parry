use crate::math::{Isometry, Real, Vector};
use crate::query::gjk::{self, CSOPoint, GJKResult, VoronoiSimplex};
use crate::shape::SupportMap;

use na::{self, Unit};
use num::Bounded;

/// Distance between support-mapped shapes.
pub fn distance_support_map_support_map<G1, G2>(pos12: &Isometry, g1: &G1, g2: &G2) -> Real
where
    G1: ?Sized + SupportMap,
    G2: ?Sized + SupportMap,
{
    distance_support_map_support_map_with_params(pos12, g1, g2, &mut VoronoiSimplex::new(), None)
}

/// Distance between support-mapped shapes.
///
/// This allows a more fine grained control other the underlying GJK algorigtm.
pub fn distance_support_map_support_map_with_params<G1, G2>(
    pos12: &Isometry,
    g1: &G1,
    g2: &G2,
    simplex: &mut VoronoiSimplex,
    init_dir: Option<Vector>,
) -> Real
where
    G1: ?Sized + SupportMap,
    G2: ?Sized + SupportMap,
{
    // TODO: or m2.translation - m1.translation ?
    let dir = init_dir.unwrap_or_else(|| -pos12.translation.vector);

    if let Some(dir) = Unit::try_new(dir, crate::math::DEFAULT_EPSILON) {
        simplex.reset(CSOPoint::from_shapes(pos12, g1, g2, &dir));
    } else {
        simplex.reset(CSOPoint::from_shapes(pos12, g1, g2, &Vector::x_axis()));
    }

    match gjk::closest_points(pos12, g1, g2, Real::max_value(), true, simplex) {
        GJKResult::Intersection => 0.0,
        GJKResult::ClosestPoints(p1, p2, _) => na::distance(&p1, &p2),
        GJKResult::Proximity(_) => unreachable!(),
        GJKResult::NoIntersection(_) => 0.0, // TODO: GJK did not converge.
    }
}
