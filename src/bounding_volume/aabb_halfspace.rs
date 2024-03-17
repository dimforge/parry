use crate::bounding_volume::Aabb;
use crate::math::*;
use crate::shape::HalfSpace;

impl HalfSpace {
    /// Computes the world-space [`Aabb`] of this half-space.
    #[inline]
    pub fn aabb(&self, _pos: &Isometry) -> Aabb {
        self.local_aabb()
    }

    /// Computes the local-space [`Aabb`] of this half-space.
    #[inline]
    pub fn local_aabb(&self) -> Aabb {
        // We divide by 2.0  so that we can still make some operations with it (like loosening)
        // without breaking the box.
        let max = Point::from(Vector::repeat(Real::MAX)) * 0.5;
        Aabb::new(-max, max)
    }
}
