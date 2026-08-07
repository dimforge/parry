use crate::bounding_volume::Aabb;
#[cfg(not(feature = "std"))]
use crate::math::ComplexField;
use crate::math::{Real, Vector};
use crate::query::{PointProjection, PointQuery, PointQueryWithLocation};
use crate::shape::{FeatureId, HeightField, TrianglePointLocation}; // For sqrt.

impl PointQuery for HeightField {
    fn project_local_point_with_max_dist(
        &self,
        pt: Vector,
        solid: bool,
        max_dist: Real,
    ) -> Option<PointProjection> {
        let aabb = Aabb::new(pt - Vector::splat(max_dist), pt + Vector::splat(max_dist));
        let mut sq_smallest_dist = Real::MAX;
        let mut best_proj = None;

        self.map_elements_in_local_aabb(&aabb, &mut |_, triangle| {
            let proj = triangle.project_local_point(pt, solid);
            let sq_dist = (pt - proj.point).length_squared();

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
    fn project_local_point(&self, point: Vector, solid: bool) -> PointProjection {
        // Grow a search neighborhood around `point` geometrically instead of iterating on
        // every element. A projection found at `dist <= max_dist` is the global closest
        // one: any closer element would intersect the ball of radius `dist` around
        // `point`, which the searched AABB contains.
        let root_aabb = self.root_aabb();
        let extents = root_aabb.extents();

        // Distance beyond which the search AABB covers every cell, making the search
        // below exhaustive.
        let max_search_dist = (point - root_aabb.center()).length() + extents.length();

        #[cfg(feature = "dim2")]
        let cell_size = self.cell_width();
        #[cfg(feature = "dim3")]
        let cell_size = self.cell_width().max(self.cell_height());

        // Initial guess: distance to the root AABB plus one cell, so the first search
        // usually visits only a few cells.
        let dist_to_aabb = (point - point.clamp(root_aabb.mins, root_aabb.maxs)).length();
        let mut search_dist = (dist_to_aabb + cell_size).max(max_search_dist * 1.0e-4);

        // TODO: the search AABB only grows, so each iteration re-tests the elements the
        //       previous ones already tested. Visit only the ring added by each step.
        while search_dist < max_search_dist {
            if let Some(proj) = self.project_local_point_with_max_dist(point, solid, search_dist) {
                return proj;
            }

            search_dist *= 4.0;
        }

        self.project_local_point_with_max_dist(point, solid, max_search_dist)
            .unwrap_or_else(|| PointProjection::new(false, point))
    }

    #[inline]
    fn project_local_point_and_get_feature(&self, point: Vector) -> (PointProjection, FeatureId) {
        // TODO: compute the feature properly.
        (self.project_local_point(point, false), FeatureId::Unknown)
    }

    // TODO: implement distance_to_point too?

    #[inline]
    fn contains_local_point(&self, _point: Vector) -> bool {
        false
    }
}

impl PointQueryWithLocation for HeightField {
    type Location = (usize, TrianglePointLocation);

    #[inline]
    fn project_local_point_and_get_location(
        &self,
        _point: Vector,
        _: bool,
    ) -> (PointProjection, Self::Location) {
        unimplemented!()
    }
}
