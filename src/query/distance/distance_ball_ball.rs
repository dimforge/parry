use crate::math::{Point, Real};
use crate::shape::Ball;
use na::{self, ComplexField};

/// Distance between balls.
#[inline]
pub fn distance_ball_ball(b1: &Ball, center2: &Point, b2: &Ball) -> Real {
    let r1 = b1.radius;
    let r2 = b2.radius;
    let distance_squared = center2.coords.norm_squared();
    let sum_radius = r1 + r2;

    if distance_squared <= sum_radius * sum_radius {
        0.0
    } else {
        ComplexField::sqrt(distance_squared) - sum_radius
    }
}
