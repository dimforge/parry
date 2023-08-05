use crate::bounding_volume::Aabb;
use crate::math::{Isometry, Point, Real};
use crate::num::Bounded;
use crate::shape::HalfSpace;
use na;

impl HalfSpace {
    /// Computes the world-space [`Aabb`] of this half-space.
    #[inline]
    pub fn aabb(&self, _pos: &Isometry<Real>) -> Aabb {
        self.local_aabb()
    }

    /// Computes the local-space [`Aabb`] of this half-space.
    #[inline]
    pub fn local_aabb(&self) -> Aabb {
        // We divide by 2.0  so that we can still make some operations with it (like loosening)
        // without breaking the box.
        let max = Point::max_value() * na::convert::<f64, Real>(0.5f64);
        Aabb::new(-max, max)
    }
}
