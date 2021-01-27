use crate::bounding_volume;
use crate::bounding_volume::BoundingSphere;
use crate::math::{Isometry, Real};
use crate::shape::ConvexPolyhedron;

impl ConvexPolyhedron {
    /// Computes the world-space bounding sphere of this convex polyhedron, transformed by `pos`.
    #[inline]
    pub fn bounding_sphere(&self, pos: &Isometry<Real>) -> BoundingSphere {
        let bv: BoundingSphere = self.local_bounding_sphere();
        bv.transform_by(pos)
    }

    /// Computes the local-space bounding sphere of this convex polyhedron.
    #[inline]
    pub fn local_bounding_sphere(&self) -> BoundingSphere {
        let (center, radius) = bounding_volume::details::point_cloud_bounding_sphere(self.points());

        BoundingSphere::new(center, radius)
    }
}
