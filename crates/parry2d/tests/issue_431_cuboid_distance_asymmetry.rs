// Regression test for https://github.com/dimforge/parry/issues/431
//
// `query::distance` between two cuboids is dispatched to a SAT-based special
// case whose face-vertex branch projected an unclamped support corner,
// overestimating the distance in face-face configurations and making the result
// depend on the arguments order. That was fixed by "fix some ambiguities in
// cuboid-cuboid SAT" (#436); these tests pin the behavior down.

use parry2d::math::{Pose, Real, Vector};
use parry2d::query::{self, ClosestPoints};
use parry2d::shape::Cuboid;

fn check_symmetric_and_exact(p1: &Pose, c1: &Cuboid, p2: &Pose, c2: &Cuboid, expected: Real) {
    let d12 = query::distance(p1, c1, p2, c2).unwrap();
    let d21 = query::distance(p2, c2, p1, c1).unwrap();

    // Cross-check against the exact GJK closest points.
    let gjk_dist = match query::closest_points(p1, c1, p2, c2, Real::MAX).unwrap() {
        ClosestPoints::WithinMargin(a, b) => (a - b).length(),
        _ => 0.0,
    };

    assert!(
        (d12 - d21).abs() < 1.0e-5,
        "asymmetric distance: {d12} vs {d21}"
    );
    assert!(
        (d12 - gjk_dist).abs() < 1.0e-5,
        "distance {d12} disagrees with GJK closest points {gjk_dist}"
    );
    assert!(
        (d12 - expected).abs() < 1.0e-4,
        "distance {d12} != expected {expected}"
    );
}

// The exact repro from issue #431: face-face configuration where the support
// corner of the taller cuboid does not project inside the other cuboid's face.
#[test]
fn issue_431_repro() {
    let c1 = Cuboid::new(Vector::new(1.0, 1.0));
    let c2 = Cuboid::new(Vector::new(1.0, 2.0));
    let p1 = Pose::identity();
    let p2 = Pose::new(Vector::new(-5.573167, 0.0), 0.0);

    check_symmetric_and_exact(&p1, &c1, &p2, &c2, 3.5731668);
}

// Separated face-face configuration (equal heights).
#[test]
fn face_face_separated() {
    let c1 = Cuboid::new(Vector::new(1.0, 1.0));
    let c2 = Cuboid::new(Vector::new(2.0, 1.0));
    let p1 = Pose::identity();
    let p2 = Pose::new(Vector::new(7.0, 0.5), 0.0);

    check_symmetric_and_exact(&p1, &c1, &p2, &c2, 4.0);
}

// Corner-corner configuration: closest features are two vertices.
#[test]
fn corner_corner_separated() {
    let c1 = Cuboid::new(Vector::new(1.0, 1.0));
    let c2 = Cuboid::new(Vector::new(1.0, 1.0));
    let p1 = Pose::identity();
    let p2 = Pose::new(Vector::new(5.0, 4.0), 0.0);

    let expected = (3.0f32 * 3.0 + 2.0 * 2.0).sqrt();
    check_symmetric_and_exact(&p1, &c1, &p2, &c2, expected);
}

// Rotated cuboid: vertex-face configuration.
#[test]
fn vertex_face_rotated() {
    let c1 = Cuboid::new(Vector::new(1.0, 1.0));
    let c2 = Cuboid::new(Vector::new(1.0, 1.0));
    let p1 = Pose::identity();
    // c2 rotated by 45 degrees, its bottom vertex facing c1's top face.
    let p2 = Pose::new(Vector::new(0.0, 4.0), core::f32::consts::FRAC_PI_4);

    let expected = 3.0 - core::f32::consts::SQRT_2;
    check_symmetric_and_exact(&p1, &c1, &p2, &c2, expected);
}

// Touching and overlapping cuboids must report a zero distance in both orders.
#[test]
fn touching_and_overlapping() {
    let c1 = Cuboid::new(Vector::new(1.0, 1.0));
    let c2 = Cuboid::new(Vector::new(1.0, 2.0));
    let p1 = Pose::identity();

    for x in [2.0, 1.5] {
        let p2 = Pose::new(Vector::new(x, 0.0), 0.0);
        let d12 = query::distance(&p1, &c1, &p2, &c2).unwrap();
        let d21 = query::distance(&p2, &c2, &p1, &c1).unwrap();
        assert!(d12.abs() < 1.0e-6, "expected zero distance, got {d12}");
        assert!(d21.abs() < 1.0e-6, "expected zero distance, got {d21}");
    }
}
