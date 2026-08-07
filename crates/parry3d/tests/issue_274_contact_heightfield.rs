// Regression test for https://github.com/dimforge/parry/issues/274
//
// `query::contact` with a heightfield used to return `Err(Unsupported)`, `HeightField`
// implementing neither `SupportMap` nor `CompositeShape`. It now iterates the elements
// intersecting the other shape's Aabb and keeps the deepest contact.

use parry3d::math::{Pose, Vector};
use parry3d::query::contact;
use parry3d::shape::{Ball, Capsule, HeightField};
use parry3d::utils::Array2;

const EPS: f32 = 1.0e-4;

fn flat_heightfield() -> HeightField {
    // A flat heightfield at y = 0 spanning x, z in [-5, 5].
    HeightField::new(
        Array2::new(3, 3, vec![0.0; 9]),
        Vector::new(10.0, 1.0, 10.0),
    )
}

#[test]
fn capsule_above_heightfield() {
    let hf = flat_heightfield();
    let capsule = Capsule::new_y(0.5, 0.2);
    let pos_hf = Pose::identity();
    // The capsule's lowest point is at y = 2.0 - 0.7 = 1.3.
    let pos_capsule = Pose::from_translation(Vector::new(0.0, 2.0, 0.0));

    // Out of prediction range: no contact, but not Unsupported.
    let c = contact(&pos_hf, &hf, &pos_capsule, &capsule, 1.0).unwrap();
    assert!(c.is_none(), "expected no contact, got {c:?}");

    // Within prediction range: a contact with the correct separation and normal.
    let c = contact(&pos_hf, &hf, &pos_capsule, &capsule, 2.0)
        .unwrap()
        .unwrap();
    assert!(
        (c.dist - 1.3).abs() < EPS,
        "expected dist ~ 1.3, got {}",
        c.dist
    );
    assert!(
        (c.normal1 - Vector::Y).length() < EPS,
        "normal1 should be +Y, got {:?}",
        c.normal1
    );
}

#[test]
fn capsule_touching_heightfield() {
    let hf = flat_heightfield();
    let capsule = Capsule::new_y(0.5, 0.2);
    let pos_hf = Pose::identity();
    // NOTE: touch the interior of a heightfield triangle: on shared triangle
    // vertices/edges the single-contact normal is ambiguous (the query sees
    // zero-thickness triangles, not the whole surface).
    let pos_capsule = Pose::from_translation(Vector::new(1.0, 0.7, -2.0));

    let c = contact(&pos_hf, &hf, &pos_capsule, &capsule, 0.1)
        .unwrap()
        .unwrap();
    assert!(c.dist.abs() < EPS, "expected touching dist, got {}", c.dist);
    assert!(
        (c.normal1 - Vector::Y).length() < EPS,
        "normal1 should be +Y, got {:?}",
        c.normal1
    );
}

#[test]
fn capsule_penetrating_heightfield_both_orders() {
    let hf = flat_heightfield();
    let capsule = Capsule::new_y(0.5, 0.2);
    let pos_hf = Pose::identity();
    // NOTE: (1, -2) is strictly inside one heightfield triangle (not on the cell
    // diagonal), and the penetration (0.1) is smaller than the capsule radius so
    // the bottom sphere's center stays strictly above the triangles' plane. Both
    // matter: on such degenerate configurations the (pre-existing) single-contact
    // query against zero-thickness triangles has an ambiguous normal.
    let pos_capsule = Pose::from_translation(Vector::new(1.0, 0.6, -2.0));

    // NOTE: EPA linearizes the capsule's round surface, hence the looser tolerance.
    const PEN_EPS: f32 = 5.0e-3;

    let c = contact(&pos_hf, &hf, &pos_capsule, &capsule, 0.0)
        .unwrap()
        .unwrap();
    assert!(
        (c.dist + 0.1).abs() < PEN_EPS,
        "expected dist ~ -0.1, got {}",
        c.dist
    );
    assert!(
        (c.normal1 - Vector::Y).length() < PEN_EPS,
        "normal1 should be +Y, got {:?}",
        c.normal1
    );

    // Flipped arguments order: the contact must be flipped as well.
    let c = contact(&pos_capsule, &capsule, &pos_hf, &hf, 0.0)
        .unwrap()
        .unwrap();
    assert!(
        (c.dist + 0.1).abs() < PEN_EPS,
        "expected dist ~ -0.1, got {}",
        c.dist
    );
    assert!(
        (c.normal1 + Vector::Y).length() < PEN_EPS,
        "normal1 should be -Y, got {:?}",
        c.normal1
    );
}

#[test]
fn ball_vs_heightfield() {
    let hf = flat_heightfield();
    let ball = Ball::new(0.5);
    let pos_hf = Pose::identity();
    let pos_ball = Pose::from_translation(Vector::new(1.0, 0.3, -2.0));

    let c = contact(&pos_hf, &hf, &pos_ball, &ball, 0.0)
        .unwrap()
        .unwrap();
    assert!(
        (c.dist + 0.2).abs() < EPS,
        "expected dist ~ -0.2, got {}",
        c.dist
    );
    assert!((c.normal1 - Vector::Y).length() < EPS);
    // The witness point on the heightfield is on its surface, under the ball.
    assert!((c.point1 - Vector::new(1.0, 0.0, -2.0)).length() < EPS);
}
