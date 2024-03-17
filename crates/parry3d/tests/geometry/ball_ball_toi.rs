// Issue #35

use parry3d::math::{GlamVectorOps, Isometry, Real, Vector};
use parry3d::query;
use parry3d::shape::Ball;

#[test]
fn test_ball_ball_toi() {
    let b = Ball::new(0.5);
    let m1 = Isometry::identity();
    let m2 = Isometry::translation(0.0, 10.0, 0.0);
    let vel1 = Vector::new(0.0, 10.0, 0.0);
    let vel2 = Vector::zeros();

    let cast = query::time_of_impact(&m1, &vel1, &b, &m2, &vel2, &b, Real::MAX, true).unwrap();

    assert_eq!(cast.unwrap().toi, 0.9);
}
