use na::{Point3, Vector3};
use parry3d::query::IntersectResult;
use parry3d::shape::TriMesh;

fn build_diamond() -> TriMesh {
    // Two tetrahedrons sharing a face
    let points = vec![
        Point3::new(0.0, 2.0, 0.0),
        Point3::new(-2.0, -1.0, 0.0),
        Point3::new(0.0, 0.0, 2.0),
        Point3::new(2.0, -1.0, 0.0),
        Point3::new(0.0, 0.0, -2.0),
    ];

    let indices = vec![
        [0u32, 1, 2],
        [0, 2, 3],
        [1, 2, 3],
        [0, 1, 4],
        [0, 4, 3],
        [1, 4, 3],
    ];

    TriMesh::new(points, indices)
}

#[test]
fn trimesh_plane_edge_intersection() {
    let mesh = build_diamond();

    let result = mesh.intersection_with_local_plane(&Vector3::ith_axis(2), 0.5, std::f32::EPSILON);

    assert!(matches!(result, IntersectResult::Intersect(_)));

    if let IntersectResult::Intersect(lines) = result {
        assert_eq!(lines.len(), 1);

        // Need to check points individually since order is not garunteed
        let vertices = lines[0].vertices();
        assert_eq!(vertices.len(), 3);
        assert!(vertices.contains(&Point3::new(-1.5, -0.75, 0.5)));
        assert!(vertices.contains(&Point3::new(1.5, -0.75, 0.5)));
        assert!(vertices.contains(&Point3::new(0.0, 1.5, 0.5)));
    }
}

#[test]
fn trimesh_plane_vertex_intersection() {
    let mesh = build_diamond();

    let result = mesh.intersection_with_local_plane(&Vector3::ith_axis(2), 0.0, std::f32::EPSILON);

    assert!(matches!(result, IntersectResult::Intersect(_)));

    if let IntersectResult::Intersect(lines) = result {
        assert_eq!(lines.len(), 1);

        // Need to check points individually since order is not garunteed
        let vertices = lines[0].vertices();
        assert_eq!(vertices.len(), 3);
        assert!(vertices.contains(&Point3::new(-2.0, -1.0, 0.0)));
        assert!(vertices.contains(&Point3::new(2.0, -1.0, 0.0)));
        assert!(vertices.contains(&Point3::new(0.0, 2.0, 0.0)));
    }
}

#[test]
fn trimesh_plane_mixed_intersection() {
    let mesh = build_diamond();

    let result = mesh.intersection_with_local_plane(&Vector3::ith_axis(0), 0.0, std::f32::EPSILON);

    assert!(matches!(result, IntersectResult::Intersect(_)));

    if let IntersectResult::Intersect(lines) = result {
        assert_eq!(lines.len(), 1);

        // Need to check points individually since order is not garunteed
        let vertices = lines[0].vertices();
        assert_eq!(vertices.len(), 4);
        assert!(vertices.contains(&Point3::new(0.0, 2.0, 0.0)));
        assert!(vertices.contains(&Point3::new(0.0, 0.0, 2.0)));
        assert!(vertices.contains(&Point3::new(0.0, -1.0, 0.0)));
        assert!(vertices.contains(&Point3::new(0.0, 0.0, -2.0)));
    }
}

#[test]
fn trimesh_plane_above() {
    let mesh = build_diamond();

    let result = mesh.intersection_with_local_plane(&Vector3::ith_axis(2), -5.0, std::f32::EPSILON);

    assert!(matches!(result, IntersectResult::Positive));
}

#[test]
fn trimesh_plane_below() {
    let mesh = build_diamond();

    let result = mesh.intersection_with_local_plane(&Vector3::ith_axis(2), 5.0, std::f32::EPSILON);

    assert!(matches!(result, IntersectResult::Negative));
}
