// Issue #35

use na::{self, Isometry3, Vector3};
use parry3d::math::Real;
use parry3d::query;
use parry3d::shape::Ball;

#[test]
fn test_ball_ball_toi() {
    let b = Ball::new(0.5);
    let m1 = Isometry3::identity();
    let m2 = Isometry3::translation(0.0, 10.0, 0.0);
    let vel1 = Vector3::new(0.0, 10.0, 0.0);
    let vel2 = Vector3::zeros();

    let cast = query::cast_shapes(&m1, &vel1, &b, &m2, &vel2, &b, Real::MAX, true).unwrap();

    assert_eq!(cast.unwrap().time_of_impact, 0.9);
}
