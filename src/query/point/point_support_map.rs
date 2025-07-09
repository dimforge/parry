use na::Unit;

use crate::math::{Isometry, Point, Real, Vector};
#[cfg(feature = "std")]
use crate::query::epa::EPA;
use crate::query::gjk::{self, CSOPoint, ConstantOrigin, GjkOptions, VoronoiSimplex};
#[cfg(feature = "dim2")]
use crate::query::point::point_query::QueryOptions;
#[cfg(feature = "dim3")]
use crate::query::point::point_query::QueryOptions;
use crate::query::{PointProjection, PointQuery};
#[cfg(feature = "dim2")]
#[cfg(feature = "std")]
use crate::shape::ConvexPolygon;
#[cfg(feature = "dim3")]
#[cfg(feature = "std")]
use crate::shape::ConvexPolyhedron;
use crate::shape::{FeatureId, SupportMap};

/// Projects a point on a shape using the GJK algorithm.
pub fn local_point_projection_on_support_map<G>(
    shape: &G,
    simplex: &mut VoronoiSimplex,
    point: &Point<Real>,
    solid: bool,
    gjk_options: &GjkOptions,
) -> PointProjection
where
    G: SupportMap,
{
    let m = Isometry::new(-point.coords, na::zero());
    let m_inv = Isometry::new(point.coords, na::zero());
    let dir = Unit::try_new(-m.translation.vector, crate::math::DEFAULT_EPSILON)
        .unwrap_or(Vector::x_axis());
    let support_point = CSOPoint::from_shapes(&m_inv, shape, &ConstantOrigin, &dir);

    simplex.reset(support_point);

    if let Some(proj) = gjk::project_origin(&m, shape, simplex, gjk_options) {
        PointProjection::new(false, proj)
    } else if solid {
        PointProjection::new(true, *point)
    } else {
        let mut epa = EPA::new();
        if let Some(pt) = epa.project_origin(&m, shape, simplex) {
            return PointProjection::new(true, pt);
        } else {
            // return match minkowski_sampling::project_origin(&m, shape, simplex) {
            //     Some(p) => PointProjection::new(true, p + point.coords),
            //     None => PointProjection::new(true, *point),
            // };

            //// All failed.
            PointProjection::new(true, *point)
        }
    }
}

#[cfg(feature = "dim3")]
impl PointQuery for ConvexPolyhedron {
    #[inline]
    fn project_local_point(
        &self,
        point: &Point<Real>,
        solid: bool,
        options: &dyn QueryOptions,
    ) -> PointProjection {
        let Some(options) = options.as_any().downcast_ref() else {
            log::warn!("Incorrect option passed to project_local_point: using default options.");
            return local_point_projection_on_support_map(
                self,
                &mut VoronoiSimplex::new(),
                point,
                solid,
                &GjkOptions::default(),
            );
        };
        local_point_projection_on_support_map(
            self,
            &mut VoronoiSimplex::new(),
            point,
            solid,
            options,
        )
    }

    #[inline]
    fn project_local_point_and_get_feature(
        &self,
        point: &Point<Real>,
        options: &dyn QueryOptions,
    ) -> (PointProjection, FeatureId) {
        let proj = self.project_local_point(point, false, options);
        let dpt = *point - proj.point;
        let local_dir = if proj.is_inside { -dpt } else { dpt };

        if let Some(local_dir) = Unit::try_new(local_dir, crate::math::DEFAULT_EPSILON) {
            let feature = self.support_feature_id_toward(&local_dir);
            (proj, feature)
        } else {
            (proj, FeatureId::Unknown)
        }
    }
}

#[cfg(feature = "dim2")]
impl PointQuery for ConvexPolygon {
    #[inline]
    fn project_local_point(
        &self,
        point: &Point<Real>,
        solid: bool,
        options: &dyn QueryOptions,
    ) -> PointProjection {
        let Some(options) = options.as_any().downcast_ref() else {
            log::warn!("Incorrect option passed to project_local_point: using default options.");
            return local_point_projection_on_support_map(
                self,
                &mut VoronoiSimplex::new(),
                point,
                solid,
                &GjkOptions::default(),
            );
        };

        return local_point_projection_on_support_map(
            self,
            &mut VoronoiSimplex::new(),
            point,
            solid,
            &options,
        );
    }

    #[inline]
    fn project_local_point_and_get_feature(
        &self,
        point: &Point<Real>,
        options: &dyn QueryOptions,
    ) -> (PointProjection, FeatureId) {
        let proj = self.project_local_point(point, false, options);
        let dpt = *point - proj.point;
        let local_dir = if proj.is_inside { -dpt } else { dpt };

        if let Some(local_dir) = Unit::try_new(local_dir, crate::math::DEFAULT_EPSILON) {
            let feature = self.support_feature_id_toward(&local_dir);
            (proj, feature)
        } else {
            (proj, FeatureId::Unknown)
        }
    }
}
