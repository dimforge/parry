use crate::bounding_volume::BoundingSphere;
use crate::math::*;
use crate::query::{PointProjection, PointQuery};
use crate::shape::{Ball, FeatureId};

impl PointQuery for BoundingSphere {
    #[inline]
    fn project_local_point(&self, pt: &Point, solid: bool) -> PointProjection {
        let centered_pt = *pt - self.center().as_vector();
        let mut proj = Ball::new(self.radius()).project_local_point(&centered_pt, solid);

        proj.point += self.center().as_vector();
        proj
    }

    #[inline]
    fn project_local_point_and_get_feature(&self, pt: &Point) -> (PointProjection, FeatureId) {
        (self.project_local_point(pt, false), FeatureId::Face(0))
    }

    #[inline]
    fn distance_to_local_point(&self, pt: &Point, solid: bool) -> Real {
        let centered_pt = *pt - self.center().as_vector();
        Ball::new(self.radius()).distance_to_local_point(&centered_pt, solid)
    }

    #[inline]
    fn contains_local_point(&self, pt: &Point) -> bool {
        let centered_pt = *pt - self.center().as_vector();
        Ball::new(self.radius()).contains_local_point(&centered_pt)
    }
}
