// Regression test for https://github.com/dimforge/parry/issues/70
// (MRP by the issue author, reformatted as a test by @ThierryBerger)
//
// Sweeping a capsule through a large cuboid produced many false negatives: the absolute
// GJK tolerance is too tight for ~100-unit support coordinates, so rounding noise made
// near-touching configurations look disjoint. The tolerance is now relative.

use parry3d::math::Pose;
use parry3d::math::Vector;
use parry3d::query::intersection_test;
use parry3d::shape::{Ball, Capsule, Cuboid, HalfSpace};

#[test]
fn capsule_cuboid_sweep_has_no_false_negatives() {
    let capsule = Capsule::new(
        Vector::new(0.0, -0.5, 0.0),
        Vector::new(0.0, 0.5, 0.0),
        0.5,
    );
    let ball = Ball::new(0.5);
    let halfspace = HalfSpace::new(Vector::new(0.0, 1.0, 0.0));

    // Upper face of the cuboid is coplanar with the outer face of the halfspace.
    let cuboid = Cuboid::new(Vector::new(50.0, 50.0, 50.0));
    let cuboid_pos = Pose::translation(0.0, -50.0, 0.0);

    let steps = 200;
    let y_max = 0.5;
    let y_min = -0.5;
    let step_size = (y_max - y_min) / steps as f32;

    let mut capsule_cuboid = 0;
    let mut capsule_halfspace = 0;
    let mut ball_cuboid = 0;
    let mut ball_halfspace = 0;

    for step in 0..steps {
        let y = y_min + step_size * step as f32;
        let test_pos = Pose::translation(0.0, y, 0.0);

        if intersection_test(&test_pos, &capsule, &Pose::IDENTITY, &halfspace).unwrap() {
            capsule_halfspace += 1;
        }
        if intersection_test(&test_pos, &capsule, &cuboid_pos, &cuboid).unwrap() {
            capsule_cuboid += 1;
        }
        if intersection_test(&test_pos, &ball, &Pose::IDENTITY, &halfspace).unwrap() {
            ball_halfspace += 1;
        }
        if intersection_test(&test_pos, &ball, &cuboid_pos, &cuboid).unwrap() {
            ball_cuboid += 1;
        }
    }

    assert_eq!(capsule_halfspace, steps);
    assert_eq!(capsule_cuboid, steps);
    assert_eq!(ball_halfspace, steps);
    assert_eq!(ball_cuboid, steps);
}
