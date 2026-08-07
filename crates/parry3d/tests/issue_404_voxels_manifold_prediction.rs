// Regression test for the ground-detection bug reported in
// https://github.com/dimforge/parry/issues/404 (comment by @bolshoytoster)
//
// `contact_manifolds_voxels_shape` did not loosen the shape's AABB by `prediction`
// before intersecting it with the voxels' AABB, so a shape hovering within the
// prediction distance above a voxels floor got no manifold at all.

use parry3d::math::{IVector, Pose, Vector};
use parry3d::query::{ContactManifold, DefaultQueryDispatcher, PersistentQueryDispatcher};
use parry3d::shape::{Ball, Capsule, Voxels};

fn voxels_floor() -> Voxels {
    // An 8x1x8 floor of unit voxels: top face at y = 1.
    let mut coords = Vec::new();
    for i in 0..8 {
        for k in 0..8 {
            coords.push(IVector::new(i, 0, k));
        }
    }
    Voxels::new(Vector::new(1.0, 1.0, 1.0), &coords)
}

fn closest_dist(manifolds: &[ContactManifold<(), ()>]) -> Option<f32> {
    manifolds
        .iter()
        .flat_map(|m| m.points.iter())
        .map(|pt| pt.dist)
        .min_by(|a, b| a.partial_cmp(b).unwrap())
}

#[test]
fn voxels_shape_manifold_within_prediction_distance() {
    let voxels = voxels_floor();

    // A capsule (the typical character-controller shape) hovering 0.05 above
    // the floor. Its lowest point is at y = 1.05.
    let capsule = Capsule::new_y(0.3, 0.2);
    let pos12 = Pose::translation(4.0, 1.55, 4.0);

    let mut manifolds: Vec<ContactManifold<(), ()>> = Vec::new();
    let mut workspace = None;
    DefaultQueryDispatcher
        .contact_manifolds(
            &pos12,
            &voxels,
            &capsule,
            0.1,
            &mut manifolds,
            &mut workspace,
        )
        .expect("the voxels/capsule pair must be supported");

    let dist = closest_dist(&manifolds)
        .expect("a speculative contact must be generated within the prediction distance");
    assert!(
        (dist - 0.05).abs() < 1.0e-3,
        "expected a contact at dist ~ 0.05, got {dist}"
    );
}

#[test]
fn voxels_shape_manifold_within_prediction_distance_direct_ball() {
    // Same check exercising `contact_manifolds_voxels_shape_shapes` directly
    // (the dispatcher would route a ball to the specialized voxels-ball
    // generator, which already loosened its AABBs correctly).
    use parry3d::query::details::contact_manifolds_voxels_shape_shapes;

    let voxels = voxels_floor();
    let ball = Ball::new(0.5);
    let pos12 = Pose::translation(4.0, 1.55, 4.0);

    let mut manifolds: Vec<ContactManifold<(), ()>> = Vec::new();
    let mut workspace = None;
    contact_manifolds_voxels_shape_shapes(
        &DefaultQueryDispatcher,
        &pos12,
        &voxels,
        &ball,
        0.1,
        &mut manifolds,
        &mut workspace,
    );

    let dist = closest_dist(&manifolds)
        .expect("a speculative contact must be generated within the prediction distance");
    assert!(
        (dist - 0.05).abs() < 1.0e-3,
        "expected a contact at dist ~ 0.05, got {dist}"
    );
}
