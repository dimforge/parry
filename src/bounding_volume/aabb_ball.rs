use crate::bounding_volume::AABB;
use crate::math::{Isometry, Point, Real, Vector};
use crate::shape::Ball;

/// Computes the Axis-Aligned Bounding Box of a ball transformed by `center`.
#[inline]
pub fn ball_aabb(center: &Point<Real>, radius: Real) -> AABB {
    AABB::new(
        *center + Vector::repeat(-radius),
        *center + Vector::repeat(radius),
    )
}

/// Computes the Axis-Aligned Bounding Box of a ball.
#[inline]
pub fn local_ball_aabb(radius: Real) -> AABB {
    let half_extents = Point::from(Vector::repeat(radius));

    AABB::new(-half_extents, half_extents)
}

impl Ball {
    /// Computes the world-space AABB of this ball transformed by `pos`.
    #[inline]
    pub fn aabb(&self, pos: &Isometry<Real>) -> AABB {
        ball_aabb(&Point::<Real>::from(pos.translation.vector), self.radius)
    }

    /// Computes the local-space AABB of this ball.
    #[inline]
    pub fn local_aabb(&self) -> AABB {
        local_ball_aabb(self.radius)
    }
}
