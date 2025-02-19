use crate::bounding_volume::BoundingSphere;
use crate::math::{Isometry, Point};
use crate::shape::Cuboid;

impl Cuboid {
    /// Computes the world-space bounding sphere of this cuboid, transformed by `pos`.
    #[inline]
    pub fn bounding_sphere(&self, pos: &Isometry) -> BoundingSphere {
        let bv: BoundingSphere = self.local_bounding_sphere();
        bv.transform_by(pos)
    }

    /// Computes the local-space bounding sphere of this cuboid.
    #[inline]
    pub fn local_bounding_sphere(&self) -> BoundingSphere {
        let radius = self.half_extents.norm();
        BoundingSphere::new(Point::origin(), radius)
    }
}
