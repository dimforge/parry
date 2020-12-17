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
    #[inline]
    pub fn aabb(&self, m: &Isometry<Real>) -> AABB {
        ball_aabb(&Point::<Real>::from(m.translation.vector), self.radius)
    }

    #[inline]
    pub fn local_aabb(&self) -> AABB {
        local_ball_aabb(self.radius)
    }
}
