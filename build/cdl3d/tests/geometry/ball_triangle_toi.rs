// Issue #123

use cdl3d::query;
use cdl3d::shape::{Ball, Triangle};
use na::{self, Isometry3, Point3, Vector3};

#[test]
fn ball_triangle_toi_infinite_loop_issue() {
    let b = Ball::new(0.375f32);
    let t = Triangle::new(
        Point3::new(0.5, -0.5, 0.0),
        Point3::new(-0.5, -0.5, 0.0),
        Point3::new(-0.5, 0.5, 0.0),
    );

    let m1 = Isometry3::translation(0.0, 0.0, 0.0);
    let m2 = Isometry3::translation(11.5, 5.5, 0.0);
    let vel1 = Vector3::new(0.0, 0.000000000000000000000000000000000000000006925, 0.0);
    let vel2 = Vector3::zeros();

    let cast = query::time_of_impact(
        &(m1.inverse() * m2),
        &(vel2 - vel1),
        &b,
        &t,
        std::f32::MAX,
        0.0,
    )
    .unwrap();

    println!("TOI: {:?}", cast);
    assert!(cast.is_none()); // The provided velocity is too small.
}
