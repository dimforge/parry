use crate::bounding_volume::BoundingSphere;
use crate::math::{Isometry, Point, Real};
use crate::shape::Cylinder;
use na::ComplexField;

impl Cylinder {
    /// Computes the world-space bounding sphere of this cylinder, transformed by `pos`.
    #[inline]
    pub fn bounding_sphere(&self, pos: &Isometry<Real>) -> BoundingSphere {
        let bv: BoundingSphere = self.local_bounding_sphere();
        bv.transform_by(pos)
    }

    /// Computes the local-space bounding sphere of this cylinder.
    #[inline]
    pub fn local_bounding_sphere(&self) -> BoundingSphere {
        let radius =
            ComplexField::sqrt(self.radius * self.radius + self.half_height * self.half_height);

        BoundingSphere::new(Point::origin(), radius)
    }
}
