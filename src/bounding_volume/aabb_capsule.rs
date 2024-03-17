use crate::bounding_volume::Aabb;
use crate::math::*;
use crate::shape::Capsule;

impl Capsule {
    /// The axis-aligned bounding box of this capsule.
    #[inline]
    pub fn aabb(&self, pos: &Isometry) -> Aabb {
        self.transform_by(pos).local_aabb()
    }

    /// The axis-aligned bounding box of this capsule.
    #[inline]
    pub fn local_aabb(&self) -> Aabb {
        let a = self.segment.a;
        let b = self.segment.b;
        let mins = a.as_vector().inf(&b.as_vector()) - Vector::repeat(self.radius);
        let maxs = a.as_vector().sup(&b.as_vector()) + Vector::repeat(self.radius);
        Aabb::new(mins.into(), maxs.into())
    }
}
