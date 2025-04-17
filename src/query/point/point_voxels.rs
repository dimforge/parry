use crate::math::{Point, Real, Vector};
use crate::query::{PointProjection, PointQuery};
use crate::shape::{Cuboid, FeatureId, VoxelType, Voxels};

impl PointQuery for Voxels {
    #[inline]
    fn project_local_point(&self, pt: &Point<Real>, solid: bool) -> PointProjection {
        // TODO: optimize this very naive implementation.
        let base_cuboid = Cuboid::new(Vector::repeat(self.voxel_size() / 2.0));
        let mut smallest_dist = Real::MAX;
        let mut result = PointProjection::new(false, *pt);

        for (_, center, data) in self.centers() {
            if data.voxel_type() != VoxelType::Empty {
                let mut candidate = base_cuboid.project_local_point(&(pt - center.coords), solid);
                candidate.point += center.coords;

                let candidate_dist = (candidate.point - pt).norm();
                if candidate_dist < smallest_dist {
                    result = candidate;
                    smallest_dist = candidate_dist;
                }
            }
        }

        result
    }

    #[inline]
    fn project_local_point_and_get_feature(
        &self,
        pt: &Point<Real>,
    ) -> (PointProjection, FeatureId) {
        // TODO: get the actual feature.
        (self.project_local_point(pt, false), FeatureId::Unknown)
    }
}
