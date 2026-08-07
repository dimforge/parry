// Regression test for https://github.com/dimforge/parry/issues/215
//
// `TriMesh::intersection_with_local_plane` used to hang forever on a flat two-triangle
// quad when the plane crossed the shared edge, the polyline extraction loop having no
// termination guard. This checks the issue's repro terminates with a sane polyline.

use parry3d::math::Vector;
use parry3d::query::IntersectResult;
use parry3d::shape::TriMesh;

#[test]
fn flat_quad_plane_intersection_terminates() {
    // Exact repro from the issue.
    let points = vec![
        Vector::new(0.0, 0.0, 0.0),
        Vector::new(0.0, 0.0, 1.0),
        Vector::new(1.0, 0.0, 0.0),
        Vector::new(1.0, 0.0, 1.0),
    ];
    let indices = vec![[0, 1, 2], [1, 3, 2]];
    let trimesh = TriMesh::new(points, indices).unwrap();

    let result = trimesh.intersection_with_local_plane(Vector::X, 0.5, 0.0005);

    match result {
        IntersectResult::Intersect(polyline) => {
            // The plane x = 0.5 cuts the unit quad along a segment of length 1.
            assert!(!polyline.vertices().is_empty());
            for pt in polyline.vertices() {
                assert!((pt.x - 0.5).abs() <= 1.0e-4, "vertex {pt:?} not on plane");
                assert!(pt.y.abs() <= 1.0e-4);
                assert!((-1.0e-4..=1.0 + 1.0e-4).contains(&pt.z));
            }

            let z_min = polyline
                .vertices()
                .iter()
                .map(|p| p.z)
                .fold(f32::MAX, f32::min);
            let z_max = polyline
                .vertices()
                .iter()
                .map(|p| p.z)
                .fold(f32::MIN, f32::max);
            assert!((z_min - 0.0).abs() <= 1.0e-4);
            assert!((z_max - 1.0).abs() <= 1.0e-4);
        }
        _ => panic!("expected an Intersect result"),
    }
}
