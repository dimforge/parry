use crate::math::{Isometry, Real, Vector};
use crate::query::gjk::{self, CSOPoint, GJKResult, VoronoiSimplex};
use crate::query::ClosestPoints;
use crate::shape::SupportMap;

use na::Unit;

/// Closest points between support-mapped shapes (`Cuboid`, `ConvexHull`, etc.)
pub fn closest_points_support_map_support_map<G1: ?Sized, G2: ?Sized>(
    pos12: &Isometry<Real>,
    g1: &G1,
    g2: &G2,
    prediction: Real,
) -> ClosestPoints
where
    G1: SupportMap,
    G2: SupportMap,
{
    match closest_points_support_map_support_map_with_params(
        pos12,
        g1,
        g2,
        prediction,
        &mut VoronoiSimplex::new(),
        None,
    ) {
        GJKResult::ClosestPoints(pt1, pt2, _) => {
            ClosestPoints::WithinMargin(pt1, pos12.inverse_transform_point(&pt2))
        }
        GJKResult::NoIntersection(_) => ClosestPoints::Disjoint,
        GJKResult::Intersection => ClosestPoints::Intersecting,
        GJKResult::Proximity(_) => unreachable!(),
    }
}

/// Closest points between support-mapped shapes (`Cuboid`, `ConvexHull`, etc.)
///
/// This allows a more fine grained control other the underlying GJK algorigtm.
pub fn closest_points_support_map_support_map_with_params<G1: ?Sized, G2: ?Sized>(
    pos12: &Isometry<Real>,
    g1: &G1,
    g2: &G2,
    prediction: Real,
    simplex: &mut VoronoiSimplex,
    init_dir: Option<Vector<Real>>,
) -> GJKResult
where
    G1: SupportMap,
    G2: SupportMap,
{
    let dir = match init_dir {
        // FIXME: or pos12.translation.vector (without the minus sign) ?
        None => -pos12.translation.vector,
        Some(dir) => dir,
    };

    if let Some(dir) = Unit::try_new(dir, crate::math::DEFAULT_EPSILON) {
        simplex.reset(CSOPoint::from_shapes(pos12, g1, g2, &dir));
    } else {
        simplex.reset(CSOPoint::from_shapes(
            pos12,
            g1,
            g2,
            &Vector::<Real>::x_axis(),
        ));
    }

    gjk::closest_points(pos12, g1, g2, prediction, true, simplex)
}
