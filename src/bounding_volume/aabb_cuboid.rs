use crate::bounding_volume::AABB;
use crate::math::{Isometry, Point, Real};
use crate::shape::Cuboid;
use crate::utils::IsometryOps;

impl Cuboid {
    /// Computes the world-space AABB of this cuboid, transformed by `pos`.
    #[inline]
    pub fn aabb(&self, pos: &Isometry<Real>) -> AABB {
        let center = Point::from(pos.translation.vector);
        let ws_half_extents = pos.absolute_transform_vector(&self.half_extents);

        AABB::from_half_extents(center, ws_half_extents)
    }

    /// Computes the local-space AABB of this cuboid.
    #[inline]
    pub fn local_aabb(&self) -> AABB {
        let half_extents = Point::from(self.half_extents);

        AABB::new(-half_extents, half_extents)
    }
}
