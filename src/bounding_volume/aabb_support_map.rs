use crate::bounding_volume;
use crate::bounding_volume::AABB;
use crate::math::{Isometry, Real};
use crate::shape::Segment;
#[cfg(feature = "dim3")]
use crate::shape::{Cone, Cylinder};

#[cfg(feature = "dim3")]
impl Cone {
    #[inline]
    pub fn aabb(&self, m: &Isometry<Real>) -> AABB {
        bounding_volume::support_map_aabb(m, self)
    }

    #[inline]
    pub fn local_aabb(&self) -> AABB {
        bounding_volume::local_support_map_aabb(self)
    }
}

#[cfg(feature = "dim3")]
impl Cylinder {
    #[inline]
    pub fn aabb(&self, m: &Isometry<Real>) -> AABB {
        bounding_volume::support_map_aabb(m, self)
    }

    #[inline]
    pub fn local_aabb(&self) -> AABB {
        bounding_volume::local_support_map_aabb(self)
    }
}

impl Segment {
    #[inline]
    pub fn aabb(&self, m: &Isometry<Real>) -> AABB {
        self.transformed(m).local_aabb()
    }

    #[inline]
    pub fn local_aabb(&self) -> AABB {
        bounding_volume::local_support_map_aabb(self)
    }
}
