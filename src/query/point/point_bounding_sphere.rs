use crate::bounding_volume::BoundingSphere;
use crate::math::{Point, Real};
use crate::query::point::point_query::QueryOptions;
use crate::query::{PointProjection, PointQuery};
use crate::shape::{Ball, FeatureId};

impl PointQuery for BoundingSphere {
    #[inline]
    fn project_local_point(
        &self,
        pt: &Point<Real>,
        solid: bool,
        options: &dyn QueryOptions,
    ) -> PointProjection {
        let centered_pt = pt - self.center().coords;
        let mut proj = Ball::new(self.radius()).project_local_point(&centered_pt, solid, options);

        proj.point += self.center().coords;
        proj
    }

    #[inline]
    fn project_local_point_and_get_feature(
        &self,
        pt: &Point<Real>,
        options: &dyn QueryOptions,
    ) -> (PointProjection, FeatureId) {
        (
            self.project_local_point(pt, false, options),
            FeatureId::Face(0),
        )
    }

    #[inline]
    fn distance_to_local_point(
        &self,
        pt: &Point<Real>,
        solid: bool,
        options: &dyn QueryOptions,
    ) -> Real {
        let centered_pt = pt - self.center().coords;
        Ball::new(self.radius()).distance_to_local_point(&centered_pt, solid, options)
    }

    #[inline]
    fn contains_local_point(&self, pt: &Point<Real>, options: &dyn QueryOptions) -> bool {
        let centered_pt = pt - self.center().coords;
        Ball::new(self.radius()).contains_local_point(&centered_pt, options)
    }
}
