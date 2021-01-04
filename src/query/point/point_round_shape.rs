use crate::math::{Point, Real};
use crate::query::gjk::VoronoiSimplex;
use crate::query::{PointProjection, PointQuery};
use crate::shape::{FeatureId, RoundShape, SupportMap};

// TODO: if PointQuery had a `project_point_with_normal` method, we could just
// call this and adjust the projected point accordingly.
impl<S: SupportMap> PointQuery for RoundShape<S> {
    #[inline]
    fn project_local_point(&self, point: &Point<Real>, solid: bool) -> PointProjection {
        crate::query::details::local_point_projection_on_support_map(
            self,
            &mut VoronoiSimplex::new(),
            point,
            solid,
        )
    }

    #[inline]
    fn project_local_point_and_get_feature(
        &self,
        point: &Point<Real>,
    ) -> (PointProjection, FeatureId) {
        (self.project_local_point(point, false), FeatureId::Unknown)
    }
}
