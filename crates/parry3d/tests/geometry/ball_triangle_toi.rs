// Issue #123

use parry3d::math::{Pose, Vector};
use parry3d::query::{self, ShapeCastOptions};
use parry3d::shape::{Ball, Triangle};

#[test]
fn ball_triangle_toi_infinite_loop_issue() {
    let b = Ball::new(0.375f32);
    let t = Triangle::new(
        Vector::new(0.5, -0.5, 0.0),
        Vector::new(-0.5, -0.5, 0.0),
        Vector::new(-0.5, 0.5, 0.0),
    );

    let m1 = Pose::translation(0.0, 0.0, 0.0);
    let m2 = Pose::translation(11.5, 5.5, 0.0);
    let vel1 = Vector::new(0.0, 0.000000000000000000000000000000000000000006925, 0.0);
    let vel2 = Vector::ZERO;

    let cast =
        query::cast_shapes(&m1, vel1, &b, &m2, vel2, &t, ShapeCastOptions::default()).unwrap();

    println!("ShapeCastHit: {:?}", cast);
    assert!(cast.is_none()); // The provided velocity is too small.
}
