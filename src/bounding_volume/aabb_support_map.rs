use crate::bounding_volume;
use crate::bounding_volume::AABB;
use crate::math::{Isometry, Real};
use crate::shape::Segment;
#[cfg(feature = "dim3")]
use crate::shape::{Cone, Cylinder};

#[cfg(feature = "dim3")]
impl Cone {
    /// Computes the world-space AABB of this cone, transformed by `pos`.
    #[inline]
    pub fn aabb(&self, pos: &Isometry<Real>) -> AABB {
        bounding_volume::details::support_map_aabb(pos, self)
    }

    /// Computes the local-space AABB of this cone.
    #[inline]
    pub fn local_aabb(&self) -> AABB {
        bounding_volume::details::local_support_map_aabb(self)
    }
}

#[cfg(feature = "dim3")]
impl Cylinder {
    /// Computes the world-space AABB of this cylinder, transformed by `pos`.
    #[inline]
    pub fn aabb(&self, pos: &Isometry<Real>) -> AABB {
        bounding_volume::details::support_map_aabb(pos, self)
    }

    /// Computes the local-space AABB of this cylinder.
    #[inline]
    pub fn local_aabb(&self) -> AABB {
        bounding_volume::details::local_support_map_aabb(self)
    }
}

impl Segment {
    /// Computes the world-space AABB of this segment, transformed by `pos`.
    #[inline]
    pub fn aabb(&self, pos: &Isometry<Real>) -> AABB {
        self.transformed(pos).local_aabb()
    }

    /// Computes the local-space AABB of this segment.
    #[inline]
    pub fn local_aabb(&self) -> AABB {
        bounding_volume::details::local_support_map_aabb(self)
    }
}
