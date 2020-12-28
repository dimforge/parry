use crate::math::{Point, Real};
use crate::query::{PointProjection, PointQuery, PointQueryWithLocation};
use crate::shape::{FeatureId, HeightField, TrianglePointLocation};
use na;
use num::Bounded;

impl PointQuery for HeightField {
    #[inline]
    fn project_local_point(&self, point: &Point<Real>, _: bool) -> PointProjection {
        let mut smallest_dist = Real::max_value();
        let mut best_proj = PointProjection::new(false, *point);

        #[cfg(feature = "dim2")]
        let iter = self.segments();
        #[cfg(feature = "dim3")]
        let iter = self.triangles();
        for elt in iter {
            let proj = elt.project_local_point(point, false);
            let dist = na::distance_squared(point, &proj.point);

            if dist < smallest_dist {
                smallest_dist = dist;
                best_proj = proj;
            }
        }

        best_proj
    }

    #[inline]
    fn project_local_point_and_get_feature(
        &self,
        point: &Point<Real>,
    ) -> (PointProjection, FeatureId) {
        // FIXME: compute the feature properly.
        (self.project_local_point(point, false), FeatureId::Unknown)
    }

    // FIXME: implement distance_to_point too?

    #[inline]
    fn contains_local_point(&self, _point: &Point<Real>) -> bool {
        false
    }
}

impl PointQueryWithLocation for HeightField {
    type Location = (usize, TrianglePointLocation);

    #[inline]
    fn project_local_point_and_get_location(
        &self,
        _point: &Point<Real>,
        _: bool,
    ) -> (PointProjection, Self::Location) {
        unimplemented!()
    }
}
