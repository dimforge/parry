use crate::bounding_volume::BoundingSphere;
use crate::math::{Isometry, Real};
use crate::shape::Capsule;

impl Capsule {
    #[inline]
    pub fn bounding_sphere(&self, m: &Isometry<Real>) -> BoundingSphere {
        self.local_bounding_sphere().transform_by(m)
    }

    #[inline]
    pub fn local_bounding_sphere(&self) -> BoundingSphere {
        let radius = self.radius + self.half_height();
        BoundingSphere::new(self.center(), radius)
    }
}
