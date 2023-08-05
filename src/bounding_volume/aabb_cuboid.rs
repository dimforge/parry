use crate::bounding_volume::Aabb;
use crate::math::{Isometry, Point, Real};
use crate::shape::Cuboid;
use crate::utils::IsometryOps;

impl Cuboid {
    /// Computes the world-space [`Aabb`] of this cuboid, transformed by `pos`.
    #[inline]
    pub fn aabb(&self, pos: &Isometry<Real>) -> Aabb {
        let center = Point::from(pos.translation.vector);
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
