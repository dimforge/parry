use crate::bounding_volume::Aabb;
use crate::math::*;
use crate::shape::Cuboid;

impl Cuboid {
    /// Computes the world-space [`Aabb`] of this cuboid, transformed by `pos`.
    #[inline]
    pub fn aabb(&self, pos: &Isometry) -> Aabb {
        let center = Point::from(pos.translation);
        let ws_half_extents = pos.absolute_transform_vector(&self.half_extents);

        Aabb::from_half_extents(center, ws_half_extents)
    }

    /// Computes the local-space [`Aabb`] of this cuboid.
    #[inline]
    pub fn local_aabb(&self) -> Aabb {
        let half_extents = Point::from(self.half_extents);

        Aabb::new(-half_extents, half_extents)
    }
}
