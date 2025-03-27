use crate::bounding_volume::BoundingSphere;
use crate::math::{Isometry, Real, Translation};
use crate::shape::{Cuboid, Voxels};

impl Voxels {
    /// Computes the world-space bounding sphere of this set of voxels, transformed by `pos`.
    #[inline]
    pub fn bounding_sphere(&self, pos: &Isometry<Real>) -> BoundingSphere {
        let shift = Translation::from(self.origin + self.extents() / 2.0);
        Cuboid::new(self.extents() / 2.0).bounding_sphere(&(pos * shift))
    }

    /// Computes the local-space bounding sphere of this set of voxels.
    #[inline]
    pub fn local_bounding_sphere(&self) -> BoundingSphere {
        Cuboid::new(self.extents() / 2.0)
            .local_bounding_sphere()
            .translated(&(self.origin.coords + self.extents() / 2.0))
    }
}
