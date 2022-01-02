// Issue #194

use na::{zero, Isometry3, Point3, Vector3};
use parry3d::math::Real;
use parry3d::query;
use parry3d::shape::TriMesh;

fn build_pyramid() -> TriMesh {
    let points = vec![
        Point3::new(0.0, 1.0, 0.0),
        Point3::new(-1.0, -0.5, 0.0),
        Point3::new(0.0, -0.5, -1.0),
        Point3::new(1.0, -0.5, 0.0),
    ];

    let indices = vec![[0u32, 1, 2], [0, 2, 3], [0, 3, 1]];

    TriMesh::new(points, indices)
}

fn do_toi_test() -> Option<Real> {
    const SPEED: Real = 100000.0;

    let shape_one = build_pyramid();
    let shape_two = build_pyramid();

    let pos_one = Vector3::new(0.0, 0.0, 0.0);
    let pos_two = Vector3::new(1000.0, 0.0, 0.0);

    let transform_one = Isometry3::new(pos_one, zero());
    let transform_two = Isometry3::new(pos_two, zero());

    let vel_one = Vector3::new(SPEED, 0.0, 0.0);
    let vel_two = Vector3::new(0.0, 0.0, 0.0);

    query::time_of_impact(
        &transform_one,
        &vel_one,
        &shape_one,
        &transform_two,
        &vel_two,
        &shape_two,
        Real::MAX,
    )
    .unwrap()
    .map(|toi| toi.toi)
}

#[test]
fn trimesh_trimesh_toi() {
    let toi = do_toi_test();
    assert_eq!(toi, Some(0.00998));
}
