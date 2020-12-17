use crate::bounding_volume;
use crate::bounding_volume::AABB;
use crate::math::Matrix;
use crate::math::{Point, Real, Real, Vector};
use crate::shape::Segment;

impl Segment {
    #[inline]
    pub fn bounding_sphere(&self, m: &Isometry<Real>) -> AABB {
        // SPEED: optimize this
        bounding_volume::support_map_aabb(m, self)
    }

    #[inline]
    pub fn local_bounding_sphere(&self) -> AABB {
        bounding_volume::local_support_map_aabb(self)
    }
}
