use crate::bounding_volume;
use crate::bounding_volume::Aabb;
use crate::math::{Isometry, Real};
use crate::shape::Segment;
#[cfg(feature = "dim3")]
use crate::shape::{Cone, Cylinder};

#[cfg(feature = "dim3")]
impl Cone {
    /// Computes the world-space [`Aabb`] of this cone, transformed by `pos`.
    #[inline]
    pub fn aabb(&self, pos: &Isometry<Real>) -> Aabb {
        bounding_volume::details::support_map_aabb(pos, self)
    }

    /// Computes the local-space [`Aabb`] of this cone.
    #[inline]
    pub fn local_aabb(&self) -> Aabb {
        bounding_volume::details::local_support_map_aabb(self)
    }
}

#[cfg(feature = "dim3")]
impl Cylinder {
    /// Computes the world-space [`Aabb`] of this cylinder, transformed by `pos`.
    #[inline]
    pub fn aabb(&self, pos: &Isometry<Real>) -> Aabb {
        bounding_volume::details::support_map_aabb(pos, self)
    }

    /// Computes the local-space [`Aabb`] of this cylinder.
    #[inline]
    pub fn local_aabb(&self) -> Aabb {
        bounding_volume::details::local_support_map_aabb(self)
    }
}

impl Segment {
    /// Computes the world-space [`Aabb`] of this segment, transformed by `pos`.
    #[inline]
    pub fn aabb(&self, pos: &Isometry<Real>) -> Aabb {
        self.transformed(pos).local_aabb()
    }

    /// Computes the local-space [`Aabb`] of this segment.
    #[inline]
    pub fn local_aabb(&self) -> Aabb {
        bounding_volume::details::local_support_map_aabb(self)
    }
}
