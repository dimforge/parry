// A ball ray cast used to lose the radius once the ray started a few thousand radii away,
// so it hit the ball's center with a NaN normal. A solid cast from the center also returned
// a NaN normal.

use parry3d::math::{Real, Vector};
use parry3d::query::{Ray, RayCast};
use parry3d::shape::Ball;

#[test]
fn distant_origin_keeps_the_radius() {
    // From 1 km a 5 cm ball used to be hit at its center.
    let ball = Ball::new(0.05);
    let ray = Ray::new(Vector::new(0.0, 1000.0, 0.0), Vector::new(0.0, -1.0, 0.0));
    let hit = ball
        .cast_local_ray_and_get_normal(&ray, Real::MAX, true)
        .unwrap();
    assert!(
        (hit.time_of_impact - 999.95).abs() < 1.0e-3,
        "{}",
        hit.time_of_impact
    );
    assert!(
        (hit.normal - Vector::Y).length() < 1.0e-5,
        "{:?}",
        hit.normal
    );

    // Off the axis, the hit point and normal are those of the sphere, not of its center.
    let ball = Ball::new(0.1);
    let ray = Ray::new(Vector::new(0.07, 300.0, 0.0), Vector::new(0.0, -1.0, 0.0));
    let hit = ball
        .cast_local_ray_and_get_normal(&ray, Real::MAX, true)
        .unwrap();
    let expected = 300.0 - (0.01 as Real - 0.0049).sqrt();
    assert!(
        (hit.time_of_impact - expected).abs() < 1.0e-4,
        "{}",
        hit.time_of_impact
    );
    assert!(
        (hit.normal - Vector::new(0.7, 0.71414286, 0.0)).length() < 1.0e-3,
        "{:?}",
        hit.normal
    );
}

#[test]
fn normals_are_never_nan() {
    let ball = Ball::new(0.05);
    for (origin, dir) in [
        (Vector::ZERO, Vector::Y),
        (Vector::ZERO, Vector::ZERO),
        (Vector::new(3000.0, 0.0, 0.0), -Vector::X),
        (Vector::new(0.0, 1.0e5, 0.0), -Vector::Y),
    ] {
        for solid in [true, false] {
            if let Some(hit) =
                ball.cast_local_ray_and_get_normal(&Ray::new(origin, dir), Real::MAX, solid)
            {
                assert!(hit.normal.is_finite(), "{origin:?} {dir:?} {solid}");
            }
        }
    }
}
