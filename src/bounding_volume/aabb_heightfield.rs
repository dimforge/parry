use crate::bounding_volume::AABB;
use crate::math::{Isometry, Real};
use crate::shape::HeightField;

impl HeightField {
    /// Computes the world-space AABB of this heightfield, transformed by `pos`.
    #[inline]
    pub fn aabb(&self, pos: &Isometry<Real>) -> AABB {
        self.root_aabb().transform_by(pos)
    }

    /// Computes the local-space AABB of this heightfield.
    #[inline]
    pub fn local_aabb(&self) -> AABB {
        self.root_aabb().clone()
    }
}
