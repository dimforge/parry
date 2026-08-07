// Verification test for https://github.com/dimforge/parry/issues/26
//
// `ShapeCastOptions::target_distance` re-added the threshold-distance parameter removed
// in early versions: the cast reports a hit as soon as the shapes come within that
// distance, i.e. earlier than a plain cast.

use parry3d::math::{Pose, Vector};
use parry3d::query::{self, ShapeCastOptions};
use parry3d::shape::Ball;

#[test]
fn target_distance_makes_the_cast_hit_earlier() {
    let b1 = Ball::new(0.5);
    let b2 = Ball::new(0.5);

    let pos1 = Pose::IDENTITY;
    let pos2 = Pose::translation(3.0, 0.0, 0.0);
    let vel1 = Vector::new(1.0, 0.0, 0.0);

    let plain = query::cast_shapes(
        &pos1,
        vel1,
        &b1,
        &pos2,
        Vector::ZERO,
        &b2,
        ShapeCastOptions::default(),
    )
    .unwrap()
    .expect("the plain cast should hit");

    let with_target_distance = query::cast_shapes(
        &pos1,
        vel1,
        &b1,
        &pos2,
        Vector::ZERO,
        &b2,
        ShapeCastOptions {
            target_distance: 0.5,
            ..Default::default()
        },
    )
    .unwrap()
    .expect("the thresholded cast should hit");

    // Plain cast: the balls touch after moving 3 - (0.5 + 0.5) = 2.
    assert!((plain.time_of_impact - 2.0).abs() < 1.0e-5);
    // With a 0.5 target distance the hit is reported 0.5 units earlier.
    assert!((with_target_distance.time_of_impact - 1.5).abs() < 1.0e-5);
    assert!(with_target_distance.time_of_impact < plain.time_of_impact);
}
