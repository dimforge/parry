// Regression test for https://github.com/dimforge/parry/issues/315
//
// `contact()` between two exactly-touching cuboids used to return arbitrary results:
// GJK's first CSO support point can land exactly on the origin, leaving EPA with a
// 0-dimensional simplex answered by a hardcoded fallback. EPA now bootstraps it and
// locates the touching face.

use parry3d::math::{Pose, Vector};
use parry3d::query::contact;
use parry3d::shape::Cuboid;

const EPS: f32 = 1.0e-5;

// The exact repro from issue #315: two unit cubes touching on the plane x = 1.
#[test]
fn exactly_touching_cuboids() {
    let box1 = Cuboid::new(Vector::new(0.5, 0.5, 0.5));
    let pos1 = Pose::from_translation(Vector::new(0.5, 0.5, 0.5));
    let pos2 = Pose::from_translation(Vector::new(1.5, 0.5, 0.5));

    let c = contact(&pos1, &box1, &pos2, &box1, 0.0).unwrap().unwrap();

    assert!(c.dist.abs() < EPS, "expected touching dist, got {}", c.dist);
    assert!(
        (c.point1.x - 1.0).abs() < EPS,
        "point1 not on the touching plane: {:?}",
        c.point1
    );
    assert!(
        (c.point2.x - 1.0).abs() < EPS,
        "point2 not on the touching plane: {:?}",
        c.point2
    );
    for pt in [c.point1, c.point2] {
        assert!(
            (0.0..=1.0).contains(&pt.y) && (0.0..=1.0).contains(&pt.z),
            "witness point outside the touching faces: {pt:?}"
        );
    }
    assert!(
        (c.normal1 - Vector::X).length() < EPS,
        "normal1 should be +X, got {:?}",
        c.normal1
    );
    assert!(
        (c.normal2 + Vector::X).length() < EPS,
        "normal2 should be -X, got {:?}",
        c.normal2
    );
}

// Same configuration with the arguments flipped: the normal must flip too.
#[test]
fn exactly_touching_cuboids_flipped() {
    let box1 = Cuboid::new(Vector::new(0.5, 0.5, 0.5));
    let pos1 = Pose::from_translation(Vector::new(0.5, 0.5, 0.5));
    let pos2 = Pose::from_translation(Vector::new(1.5, 0.5, 0.5));

    let c = contact(&pos2, &box1, &pos1, &box1, 0.0).unwrap().unwrap();

    assert!(c.dist.abs() < EPS, "expected touching dist, got {}", c.dist);
    assert!((c.point1.x - 1.0).abs() < EPS);
    assert!((c.point2.x - 1.0).abs() < EPS);
    assert!(
        (c.normal1 + Vector::X).length() < EPS,
        "normal1 should be -X, got {:?}",
        c.normal1
    );
}

// Overlapping cuboids must keep returning a negative dist consistent with the
// actual penetration.
#[test]
fn overlapping_cuboids() {
    let box1 = Cuboid::new(Vector::new(0.5, 0.5, 0.5));
    let pos1 = Pose::from_translation(Vector::new(0.5, 0.5, 0.5));
    let pos2 = Pose::from_translation(Vector::new(1.2, 0.5, 0.5));

    let c = contact(&pos1, &box1, &pos2, &box1, 0.0).unwrap().unwrap();

    assert!(
        (c.dist + 0.3).abs() < 1.0e-4,
        "expected dist ~ -0.3, got {}",
        c.dist
    );
    assert!((c.normal1 - Vector::X).length() < EPS);
    // The witness points must be on the shapes' penetrating faces.
    assert!((c.point1.x - 1.0).abs() < 1.0e-4);
    assert!((c.point2.x - 0.7).abs() < 1.0e-4);
}
