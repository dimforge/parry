use crate::bounding_volume::Aabb;
use crate::math::{Isometry, Real};
use crate::shape::ConvexPolygon;

impl ConvexPolygon {
    /// Computes the world-space [`Aabb`] of this convex polygon, transformed by `pos`.
    #[inline]
    pub fn aabb(&self, pos: &Isometry<Real>) -> Aabb {
        super::details::point_cloud_aabb(pos, self.points())
    }

    /// Computes the local-space [`Aabb`] of this convex polygon.
    #[inline]
    pub fn local_aabb(&self) -> Aabb {
        super::details::local_point_cloud_aabb(self.points())
    }
}
