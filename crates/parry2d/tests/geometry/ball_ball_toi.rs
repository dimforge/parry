// Issue #35

use na::{self, Isometry2, Vector2};
use parry2d::math::Real;
use parry2d::query;
use parry2d::shape::Ball;

#[test]
fn test_ball_ball_toi() {
    let b = Ball::new(0.5);
    let m1 = Isometry2::identity();
    let m2 = Isometry2::translation(0.0, 10.0);
    let v1 = Vector2::new(0.0, 10.0);
    let v2 = Vector2::zeros();

    let cast = query::time_of_impact(&m1, &v1, &b, &m2, &v2, &b, Real::MAX).unwrap();

    assert_eq!(cast.unwrap().toi, 0.9);
}
