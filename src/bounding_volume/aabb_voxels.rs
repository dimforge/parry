use crate::bounding_volume::Aabb;
use crate::math::{Isometry, Real, Translation};
use crate::shape::{Cuboid, Voxels};

impl Voxels {
    /// Computes the world-space Aabb of this set of voxels, transformed by `pos`.
    #[inline]
    pub fn aabb(&self, pos: &Isometry<Real>) -> Aabb {
        let shift = Translation::from(self.origin + self.extents() / 2.0);
        Cuboid::new(self.extents() / 2.0).aabb(&(pos * shift))
    }

    /// Computes the local-space Aabb of this set of voxels.
    #[inline]
    pub fn local_aabb(&self) -> Aabb {
        Cuboid::new(self.extents() / 2.0)
            .local_aabb()
            .translated(&(self.origin.coords + self.extents() / 2.0))
    }
}
