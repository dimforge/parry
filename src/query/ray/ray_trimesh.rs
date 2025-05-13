use crate::math::Real;
use crate::query::{Ray, RayCast, RayIntersection};
use crate::shape::{FeatureId, TriMesh};
use crate::query::ray::{RayCompositeShapeToiAndNormalBestFirstVisitor, RayCompositeShapeToiBestFirstVisitor};

impl RayCast for TriMesh {
    #[inline]
    fn cast_local_ray(&self, ray: &Ray, max_time_of_impact: Real, solid: bool) -> Option<Real> {
        let mut visitor =
            RayCompositeShapeToiBestFirstVisitor::new(self, ray, max_time_of_impact, solid);

        self.qbvh()
            .traverse_best_first(&mut visitor)
            .map(|res| res.1 .1)
    }

    #[inline]
    fn cast_local_ray_and_get_normal(
        &self,
        ray: &Ray,
        max_time_of_impact: Real,
        solid: bool,
    ) -> Option<RayIntersection> {
        let mut visitor = RayCompositeShapeToiAndNormalBestFirstVisitor::new(
            self,
            ray,
            max_time_of_impact,
            solid,
        );

        self.qbvh()
            .traverse_best_first(&mut visitor)
            .map(|(_, (best, mut res))| {
                // We hit a backface.
                // NOTE: we need this for `TriMesh::is_backface` to work properly.
                if res.feature == FeatureId::Face(1) {
                    res.feature = FeatureId::Face(best + self.indices().len() as u32)
                } else {
                    res.feature = FeatureId::Face(best);
                }
                res
            })
    }
}