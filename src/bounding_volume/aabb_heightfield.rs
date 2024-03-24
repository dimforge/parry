use crate::bounding_volume::Aabb;
use crate::math::{Isometry, Real};
use crate::shape::{GenericHeightField, HeightFieldStorage};

impl<Storage: HeightFieldStorage> GenericHeightField<Storage> {
    /// Computes the world-space [`Aabb`] of this heightfield, transformed by `pos`.
    #[inline]
    pub fn aabb(&self, pos: &Isometry<Real>) -> Aabb {
        self.root_aabb().transform_by(pos)
    }

    /// Computes the local-space [`Aabb`] of this heightfield.
    #[inline]
    pub fn local_aabb(&self) -> Aabb {
        *self.root_aabb()
    }
}
