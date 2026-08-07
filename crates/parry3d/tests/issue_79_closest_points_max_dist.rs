// Regression test for https://github.com/dimforge/parry/issues/79
//
// `closest_points` used to panic when a `Cylinder` and a `TriMesh` were
// farther apart than `max_dist` (the old best-first visitor asserted that a
// result was always found). It must return `Ok(ClosestPoints::Disjoint)`.

use parry3d::math::{Pose, Vector};
use parry3d::query::{self, ClosestPoints};
use parry3d::shape::{Cylinder, TriMesh};

fn test_shapes() -> (Cylinder, TriMesh) {
    let cylinder = Cylinder::new(0.5, 1.0);

    // Exact mesh from the issue: one triangle at distance 0.5 below the cylinder.
    let vertices = vec![
        Vector::new(1.0, -1.0, 0.0),
        Vector::new(-1.0, -1.0, 0.0),
        Vector::new(0.0, -1.0, 1.0),
    ];
    let indices = vec![[0, 1, 2]];
    let trimesh = TriMesh::new(vertices, indices).unwrap();

    (cylinder, trimesh)
}

#[test]
fn closest_points_beyond_max_dist_is_disjoint() {
    let (cylinder, trimesh) = test_shapes();

    // The shapes are at distance 0.5: with max_dist = 0.49 this used to panic.
    let res = query::closest_points(&Pose::IDENTITY, &cylinder, &Pose::IDENTITY, &trimesh, 0.49);
    assert!(matches!(res, Ok(ClosestPoints::Disjoint)));

    // Same in the other argument order.
    let res = query::closest_points(&Pose::IDENTITY, &trimesh, &Pose::IDENTITY, &cylinder, 0.49);
    assert!(matches!(res, Ok(ClosestPoints::Disjoint)));
}

#[test]
fn closest_points_within_max_dist_is_within_margin() {
    let (cylinder, trimesh) = test_shapes();

    let res = query::closest_points(&Pose::IDENTITY, &cylinder, &Pose::IDENTITY, &trimesh, 0.6);
    match res {
        Ok(ClosestPoints::WithinMargin(p1, p2)) => {
            let dist = (p1 - p2).length();
            assert!((dist - 0.5).abs() < 1.0e-4, "distance was {dist}");
        }
        other => panic!("expected WithinMargin, got {other:?}"),
    }
}
