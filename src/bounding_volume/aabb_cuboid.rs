use crate::bounding_volume::AABB;
use crate::math::{Isometry, Point, Real};
use crate::shape::Cuboid;
use crate::utils::IsometryOps;

impl Cuboid {
    #[inline]
    pub fn aabb(&self, m: &Isometry<Real>) -> AABB {
        let center = Point::from(m.translation.vector);
        let ws_half_extents = m.absolute_transform_vector(&self.half_extents);

        AABB::from_half_extents(center, ws_half_extents)
    }

    #[inline]
    pub fn local_aabb(&self) -> AABB {
        let half_extents = Point::from(self.half_extents);

        AABB::new(-half_extents, half_extents)
    }
}
