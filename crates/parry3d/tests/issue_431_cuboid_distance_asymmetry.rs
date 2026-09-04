// Regression test for https://github.com/dimforge/parry/issues/431 (3D variant)
//
// `query::distance` between two cuboids is dispatched to a SAT-based special
// case whose face-vertex branch projected an unclamped support corner,
// overestimating the distance in face-face configurations and making the result
// depend on the arguments order. That was fixed by "fix some ambiguities in
// cuboid-cuboid SAT" (#436); these tests pin the behavior down.

use parry3d::math::{Pose, Real, Vector};
use parry3d::query::{self, ClosestPoints};
use parry3d::shape::Cuboid;

fn check_symmetric_and_exact(p1: &Pose, c1: &Cuboid, p2: &Pose, c2: &Cuboid, expected: Real) {
    let d12 = query::distance(p1, c1, p2, c2).unwrap().distance;
    let d21 = query::distance(p2, c2, p1, c1).unwrap().distance;

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

// 3D analog of the issue #431 repro: face-face configuration where the support
// corner of the larger cuboid does not project inside the other cuboid's face.
#[test]
fn face_face_overhang() {
    let c1 = Cuboid::new(Vector::new(1.0, 1.0, 1.0));
    let c2 = Cuboid::new(Vector::new(1.0, 2.0, 2.0));
    let p1 = Pose::identity();
    let p2 = Pose::from_translation(Vector::new(-5.573167, 0.0, 0.0));

    check_symmetric_and_exact(&p1, &c1, &p2, &c2, 3.5731668);
}

// Face-vertex configuration: c2 is rotated 45 degrees about Z on the diagonal,
// which aligns one of its faces with the diagonal direction; the closest pair is
// c1's corner (1, 1, z) against the interior of that face.
#[test]
fn vertex_face_diagonal() {
    let c1 = Cuboid::new(Vector::new(1.0, 1.0, 1.0));
    let c2 = Cuboid::new(Vector::new(1.0, 1.0, 1.0));
    let p1 = Pose::identity();
    let p2 = Pose::new(
        Vector::new(4.0, 4.0, 0.0),
        Vector::new(0.0, 0.0, core::f32::consts::FRAC_PI_4),
    );

    // Distance from c1's corner to c2's face plane: 3 * sqrt(2) - 1.
    let expected = 3.0 * core::f32::consts::SQRT_2 - 1.0;
    check_symmetric_and_exact(&p1, &c1, &p2, &c2, expected);
}

// Edge-edge configuration: c2 is rotated 30 degrees about Y so that one of its
// y-parallel edges faces one of c1's y-parallel edges.
#[test]
fn edge_edge_separated() {
    let c1 = Cuboid::new(Vector::new(1.0, 1.0, 1.0));
    let c2 = Cuboid::new(Vector::new(1.0, 1.0, 1.0));
    let p1 = Pose::identity();
    let angle = core::f32::consts::FRAC_PI_6;
    let p2 = Pose::new(Vector::new(5.0, 0.0, 5.0), Vector::new(0.0, angle, 0.0));

    // Closest features: c1's edge through (1, y, 1) and c2's edge through the
    // rotated image of its local (-1, y, -1) corner.
    let (s, c) = angle.sin_cos();
    let expected = ((4.0 - c - s) * (4.0 - c - s) + (4.0 - c + s) * (4.0 - c + s)).sqrt();
    check_symmetric_and_exact(&p1, &c1, &p2, &c2, expected);
}

// Corner-corner configuration.
#[test]
fn corner_corner_separated() {
    let c1 = Cuboid::new(Vector::new(1.0, 1.0, 1.0));
    let c2 = Cuboid::new(Vector::new(1.0, 1.0, 1.0));
    let p1 = Pose::identity();
    let p2 = Pose::from_translation(Vector::new(5.0, 4.0, 3.0));

    let expected = (3.0f32 * 3.0 + 2.0 * 2.0 + 1.0).sqrt();
    check_symmetric_and_exact(&p1, &c1, &p2, &c2, expected);
}

// Touching and overlapping cuboids must report a zero distance in both orders.
#[test]
fn touching_and_overlapping() {
    let c1 = Cuboid::new(Vector::new(1.0, 1.0, 1.0));
    let c2 = Cuboid::new(Vector::new(1.0, 2.0, 2.0));
    let p1 = Pose::identity();

    for x in [2.0, 1.5] {
        let p2 = Pose::from_translation(Vector::new(x, 0.0, 0.0));
        let d12 = query::distance(&p1, &c1, &p2, &c2).unwrap().distance;
        let d21 = query::distance(&p2, &c2, &p1, &c1).unwrap().distance;
        assert!(d12.abs() < 1.0e-6, "expected zero distance, got {d12}");
        assert!(d21.abs() < 1.0e-6, "expected zero distance, got {d21}");
    }
}
