use na::{self, Unit};

use crate::math::{Isometry, Real, Vector};
use crate::query::gjk::{self, CSOPoint, GJKResult, VoronoiSimplex};
use crate::shape::SupportMap;

/// Intersection test between support-mapped shapes (`Cuboid`, `ConvexHull`, etc.)
pub fn intersection_test_support_map_support_map<G1: ?Sized, G2: ?Sized>(
    pos12: &Isometry<Real>,
    g1: &G1,
    g2: &G2,
) -> bool
where
    G1: SupportMap,
    G2: SupportMap,
{
    intersection_test_support_map_support_map_with_params(
        pos12,
        g1,
        g2,
        &mut VoronoiSimplex::new(),
        None,
    )
    .0
}

/// Intersection test between support-mapped shapes (`Cuboid`, `ConvexHull`, etc.)
///
/// This allows a more fine grained control other the underlying GJK algorithm.
pub fn intersection_test_support_map_support_map_with_params<G1: ?Sized, G2: ?Sized>(
    pos12: &Isometry<Real>,
    g1: &G1,
    g2: &G2,
    simplex: &mut VoronoiSimplex,
    init_dir: Option<Unit<Vector<Real>>>,
) -> (bool, Unit<Vector<Real>>)
where
    G1: SupportMap,
    G2: SupportMap,
{
    let dir = if let Some(init_dir) = init_dir {
        init_dir
    } else if let Some(init_dir) =
        Unit::try_new(pos12.translation.vector, crate::math::DEFAULT_EPSILON)
    {
        init_dir
    } else {
        Vector::x_axis()
    };

    simplex.reset(CSOPoint::from_shapes(pos12, g1, g2, &dir));

    match gjk::closest_points(pos12, g1, g2, 0.0, false, simplex) {
        GJKResult::Intersection => (true, dir),
        GJKResult::Proximity(dir) => (false, dir),
        GJKResult::NoIntersection(dir) => (false, dir),
        GJKResult::ClosestPoints(..) => unreachable!(),
    }
}
