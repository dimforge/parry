use crate::bounding_volume::Aabb;
use crate::math::{Isometry, Real};
use crate::shape::ConvexPolyhedron;

impl ConvexPolyhedron {
    /// Computes the world-space [`Aabb`] of this convex polyhedron, transformed by `pos`.
    #[inline]
    pub fn aabb(&self, pos: &Isometry<Real>) -> Aabb {
        super::details::point_cloud_aabb_ref(pos, self.points())
    }

    /// Computes the local-space [`Aabb`] of this convex polyhedron.
    #[inline]
    pub fn local_aabb(&self) -> Aabb {
        super::details::local_point_cloud_aabb_ref(self.points())
    }
}
