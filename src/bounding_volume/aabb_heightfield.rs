use crate::bounding_volume::AABB;
use crate::math::{Isometry, Real};
use crate::shape::HeightField;

impl HeightField {
    #[inline]
    pub fn aabb(&self, m: &Isometry<Real>) -> AABB {
        self.root_aabb().transform_by(m)
    }

    #[inline]
    pub fn local_aabb(&self) -> AABB {
        self.root_aabb().clone()
    }
}
