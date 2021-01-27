use crate::bounding_volume::BoundingSphere;
use crate::math::{Isometry, Real};
use crate::shape::Capsule;

impl Capsule {
    /// Computes the world-space bounding sphere of this capsule, transformed by `pos`.
    #[inline]
    pub fn bounding_sphere(&self, pos: &Isometry<Real>) -> BoundingSphere {
        self.local_bounding_sphere().transform_by(pos)
    }

    /// Computes the world-space bounding sphere of this capsule.
    #[inline]
    pub fn local_bounding_sphere(&self) -> BoundingSphere {
        let radius = self.radius + self.half_height();
        BoundingSphere::new(self.center(), radius)
    }
}
