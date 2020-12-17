use crate::bounding_volume::aabb_utils;
use crate::bounding_volume::AABB;
use crate::math::{Isometry, Real};
use crate::shape::ConvexPolygon;

impl ConvexPolygon {
    #[inline]
    pub fn aabb(&self, m: &Isometry<Real>) -> AABB {
        aabb_utils::point_cloud_aabb(m, self.points())
    }

    #[inline]
    pub fn local_aabb(&self) -> AABB {
        aabb_utils::local_point_cloud_aabb(self.points())
    }
}
