use crate::bounding_volume::Aabb;
use crate::math::{Isometry, Real, Vector};
use crate::shape::Capsule;

impl Capsule {
    /// The axis-aligned bounding box of this capsule.
    #[inline]
    pub fn aabb(&self, pos: &Isometry<Real>) -> Aabb {
        self.transform_by(pos).local_aabb()
    }

    /// The axis-aligned bounding box of this capsule.
    #[inline]
    pub fn local_aabb(&self) -> Aabb {
        let a = self.segment.a;
        let b = self.segment.b;
        let mins = a.coords.inf(&b.coords) - Vector::repeat(self.radius);
        let maxs = a.coords.sup(&b.coords) + Vector::repeat(self.radius);
        Aabb::new(mins.into(), maxs.into())
    }
}
