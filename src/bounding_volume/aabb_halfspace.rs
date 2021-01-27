use crate::bounding_volume::AABB;
use crate::math::{Isometry, Point, Real};
use crate::num::Bounded;
use crate::shape::HalfSpace;
use na;

impl HalfSpace {
    /// Computes the world-space AABB of this half-space.
    #[inline]
    pub fn aabb(&self, _pos: &Isometry<Real>) -> AABB {
        self.local_aabb()
    }

    /// Computes the local-space AABB of this half-space.
    #[inline]
    pub fn local_aabb(&self) -> AABB {
        // We divide by 2.0  so that we can still make some operations with it (like loosening)
        // without breaking the box.
        let max = Point::max_value() * na::convert::<f64, Real>(0.5f64);
        AABB::new(-max, max)
    }
}
