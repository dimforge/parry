use crate::bounding_volume::AABB;
use crate::math::{Isometry, Real};
use crate::shape::ConvexPolyhedron;

impl ConvexPolyhedron {
    #[inline]
    pub fn aabb(&self, m: &Isometry<Real>) -> AABB {
        super::point_cloud_aabb(m, self.points())
    }

    #[inline]
    pub fn local_aabb(&self) -> AABB {
        super::local_point_cloud_aabb(self.points())
    }
}
