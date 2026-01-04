// Issue #35

use parry3d::math::{Pose, Vector};
use parry3d::query::{self, ShapeCastOptions};
use parry3d::shape::Ball;

#[test]
fn test_ball_ball_toi() {
    let b = Ball::new(0.5);
    let m1 = Pose::identity();
    let m2 = Pose::translation(0.0, 10.0, 0.0);
    let vel1 = Vector::new(0.0, 10.0, 0.0);
    let vel2 = Vector::ZERO;

    let cast =
        query::cast_shapes(&m1, vel1, &b, &m2, vel2, &b, ShapeCastOptions::default()).unwrap();

    assert_eq!(cast.unwrap().time_of_impact, 0.9);
}
