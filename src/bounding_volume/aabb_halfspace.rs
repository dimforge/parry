use crate::bounding_volume::AABB;
use crate::math::{Isometry, Point, Real};
use crate::num::Bounded;
use crate::shape::HalfSpace;
use na;

impl HalfSpace {
    #[inline]
    pub fn aabb(&self, _: &Isometry<Real>) -> AABB {
        self.local_aabb()
    }

    #[inline]
    pub fn local_aabb(&self) -> AABB {
        // We divide by 2.0  so that we can still make some operations with it (like loosening)
        // without breaking the box.
        let max = Point::max_value() * na::convert::<f64, Real>(0.5f64);
        AABB::new(-max, max)
    }
}
