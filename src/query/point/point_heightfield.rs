use crate::bounding_volume::Aabb;
use crate::math::{Point, Real, Vector};
use crate::query::{PointProjection, PointQuery, PointQueryWithLocation, QueryOptions};
use crate::shape::{FeatureId, HeightField, TrianglePointLocation};
#[cfg(not(feature = "std"))]
use na::ComplexField; // For sqrt.

impl PointQuery for HeightField {
    fn project_local_point_with_max_dist(
        &self,
        pt: &Point<Real>,
        solid: bool,
        max_dist: Real,
        options: &dyn QueryOptions,
    ) -> Option<PointProjection> {
        let aabb = Aabb::new(pt - Vector::repeat(max_dist), pt + Vector::repeat(max_dist));
        let mut sq_smallest_dist = Real::MAX;
        let mut best_proj = None;

        self.map_elements_in_local_aabb(&aabb, &mut |_, triangle| {
            let proj = triangle.project_local_point(pt, solid, options);
            let sq_dist = na::distance_squared(pt, &proj.point);

            if sq_dist < sq_smallest_dist {
                sq_smallest_dist = sq_dist;

                if sq_dist.sqrt() <= max_dist {
                    best_proj = Some(proj);
                }
            }
        });

        best_proj
    }

    #[inline]
    fn project_local_point(
        &self,
        point: &Point<Real>,
        _: bool,
        options: &dyn QueryOptions,
    ) -> PointProjection {
        let mut smallest_dist = Real::MAX;
        let mut best_proj = PointProjection::new(false, *point);

        #[cfg(feature = "dim2")]
        let iter = self.segments();
        #[cfg(feature = "dim3")]
        let iter = self.triangles();
        for elt in iter {
            let proj = elt.project_local_point(point, false, options);
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
        options: &dyn QueryOptions,
    ) -> (PointProjection, FeatureId) {
        // TODO: compute the feature properly.
        (
            self.project_local_point(point, false, options),
            FeatureId::Unknown,
        )
    }

    // TODO: implement distance_to_point too?

    #[inline]
    fn contains_local_point(&self, _point: &Point<Real>, _options: &dyn QueryOptions) -> bool {
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
        _options: &dyn QueryOptions,
    ) -> (PointProjection, Self::Location) {
        unimplemented!()
    }
}
