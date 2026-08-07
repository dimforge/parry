// Regression test for https://github.com/dimforge/parry/issues/17
//
// A cube `TriMesh` scaled down to 1e-10 used to panic on an `Option::unwrap` inside the
// old Qbvh best-first traversal, which has since been replaced (see also issue #395).

use parry3d::math::Vector;
use parry3d::query::PointQuery;
use parry3d::shape::TriMesh;

#[test]
fn tiny_trimesh_distance_to_point_does_not_panic() {
    let length_unit = 1.0e-10_f32;

    // Exact mesh from the issue.
    let vertices: Vec<_> = [
        [-0.5, -0.5, 0.5],
        [0.5, -0.5, 0.5],
        [-0.5, 0.5, 0.5],
        [0.5, 0.5, 0.5],
        [-0.5, 0.5, -0.5],
        [0.5, 0.5, -0.5],
        [-0.5, -0.5, -0.5],
        [0.5, -0.5, -0.5],
    ]
    .iter()
    .map(|p| Vector::new(p[0] * length_unit, p[1] * length_unit, p[2] * length_unit))
    .collect();
    let indices = vec![
        [3, 1, 0],
        [2, 3, 0],
        [5, 3, 2],
        [4, 5, 2],
        [7, 5, 4],
        [6, 7, 4],
        [1, 7, 6],
        [0, 1, 6],
        [5, 7, 1],
        [3, 5, 1],
        [2, 0, 6],
        [4, 2, 6],
    ];

    let trimesh = TriMesh::new(vertices, indices).unwrap();

    let queries = [
        Vector::ZERO,
        Vector::splat(0.25 * length_unit),
        Vector::splat(2.0 * length_unit),
        Vector::new(1.0, -2.0, 3.0), // Far away compared to the mesh scale.
    ];

    for pt in queries {
        let dist = trimesh.distance_to_local_point(pt, true);
        assert!(dist.is_finite());
    }
}
