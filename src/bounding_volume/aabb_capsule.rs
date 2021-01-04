use crate::bounding_volume::AABB;
use crate::math::{Isometry, Real, Vector};
use crate::shape::Capsule;

impl Capsule {
    /// The axis-aligned bounding box of this capsule.
    #[inline]
    pub fn aabb(&self, pos: &Isometry<Real>) -> AABB {
        self.transform_by(pos).local_aabb()
    }

    /// The axis-aligned bounding box of this capsule.
    #[inline]
    pub fn local_aabb(&self) -> AABB {
        let a = self.segment.a;
        let b = self.segment.b;
        let mins = a.coords.inf(&b.coords) - Vector::repeat(self.radius);
        let maxs = a.coords.sup(&b.coords) + Vector::repeat(self.radius);
        AABB::new(mins.into(), maxs.into())
    }
}
