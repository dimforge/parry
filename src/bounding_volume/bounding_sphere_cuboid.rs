use crate::bounding_volume::BoundingSphere;
use crate::math::{Isometry, Point, Real};
use crate::shape::Cuboid;

impl Cuboid {
    #[inline]
    pub fn bounding_sphere(&self, m: &Isometry<Real>) -> BoundingSphere {
        let bv: BoundingSphere = self.local_bounding_sphere();
        bv.transform_by(m)
    }

    #[inline]
    pub fn local_bounding_sphere(&self) -> BoundingSphere {
        let radius = self.half_extents.norm();
        BoundingSphere::new(Point::origin(), radius)
    }
}
