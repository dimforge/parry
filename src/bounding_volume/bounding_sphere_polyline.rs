use crate::bounding_volume::BoundingSphere;
use crate::math::{Isometry, Real};
use crate::shape::Polyline;

impl Polyline {
    /// Computes the world-space bounding sphere of this polyline, transformed by `pos`.
    #[inline]
    pub fn bounding_sphere(&self, pos: &Isometry<Real>) -> BoundingSphere {
        self.local_aabb().bounding_sphere().transform_by(pos)
    }

    /// Computes the local-space bounding sphere of this polyline.
    #[inline]
    pub fn local_bounding_sphere(&self) -> BoundingSphere {
        self.local_aabb().bounding_sphere()
    }
}
