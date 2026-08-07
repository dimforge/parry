// Regression test for https://github.com/dimforge/parry/issues/315 (2D variant)
//
// `contact()` between two exactly-touching cuboids used to return both witness points
// at the shapes' centers, fabricated by EPA's degenerate 0-dimensional simplex
// fallback. It now returns the actual witness points.

use parry2d::math::{Pose, Vector};
use parry2d::query::contact;
use parry2d::shape::Cuboid;

const EPS: f32 = 1.0e-5;

#[test]
fn exactly_touching_cuboids() {
    let box1 = Cuboid::new(Vector::new(0.5, 0.5));
    let pos1 = Pose::from_translation(Vector::new(0.5, 0.5));
    let pos2 = Pose::from_translation(Vector::new(1.5, 0.5));

    let c = contact(&pos1, &box1, &pos2, &box1, 0.0).unwrap().unwrap();

    assert!(c.dist.abs() < EPS, "expected touching dist, got {}", c.dist);
    assert!(
        (c.point1.x - 1.0).abs() < EPS,
        "point1 not on the touching line: {:?}",
        c.point1
    );
    assert!(
        (c.point2.x - 1.0).abs() < EPS,
        "point2 not on the touching line: {:?}",
        c.point2
    );
    assert!(
        (c.normal1 - Vector::X).length() < EPS,
        "normal1 should be +X, got {:?}",
        c.normal1
    );
}

#[test]
fn overlapping_cuboids() {
    let box1 = Cuboid::new(Vector::new(0.5, 0.5));
    let pos1 = Pose::from_translation(Vector::new(0.5, 0.5));
    let pos2 = Pose::from_translation(Vector::new(1.2, 0.5));

    let c = contact(&pos1, &box1, &pos2, &box1, 0.0).unwrap().unwrap();

    assert!(
        (c.dist + 0.3).abs() < 1.0e-4,
        "expected dist ~ -0.3, got {}",
        c.dist
    );
    assert!((c.normal1 - Vector::X).length() < EPS);
}
