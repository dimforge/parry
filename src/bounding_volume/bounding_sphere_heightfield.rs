use crate::bounding_volume::BoundingSphere;
use crate::math::Isometry;
use crate::shape::HeightField;

impl HeightField {
    /// Computes the world-space bounding sphere of this height-field, transformed by `pos`.
    #[inline]
    pub fn bounding_sphere(&self, pos: &Isometry) -> BoundingSphere {
        self.local_aabb().bounding_sphere().transform_by(pos)
    }

    /// Computes the local-space bounding sphere of this height-field.
    #[inline]
    pub fn local_bounding_sphere(&self) -> BoundingSphere {
        self.local_aabb().bounding_sphere()
    }
}
