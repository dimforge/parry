use parry3d::math::{Pose, Vector};
use parry3d::query::IntersectResult;
use parry3d::shape::TriMesh;

fn build_diamond(position: &Pose) -> TriMesh {
    // Two tetrahedrons sharing a face
    let points = vec![
        position * Vector::new(0.0, 2.0, 0.0),
        position * Vector::new(-2.0, -1.0, 0.0),
        position * Vector::new(0.0, 0.0, 2.0),
        position * Vector::new(2.0, -1.0, 0.0),
        position * Vector::new(0.0, 0.0, -2.0),
    ];

    let indices = vec![
        [0u32, 1, 2],
        [0, 2, 3],
        [1, 2, 3],
        [0, 1, 4],
        [0, 4, 3],
        [1, 4, 3],
    ];

    TriMesh::new(points, indices).unwrap()
}

#[test]
fn trimesh_plane_edge_intersection() {
    let mesh = build_diamond(&Pose::identity());

    let result = mesh.intersection_with_local_plane(Vector::Z, 0.5, core::f32::EPSILON);

    assert!(matches!(result, IntersectResult::Intersect(_)));

    if let IntersectResult::Intersect(line) = result {
        // Need to check points individually since order is not guaranteed
        let vertices = line.vertices();
        assert_eq!(vertices.len(), 3);
        assert!(vertices.contains(&Vector::new(-1.5, -0.75, 0.5)));
        assert!(vertices.contains(&Vector::new(1.5, -0.75, 0.5)));
        assert!(vertices.contains(&Vector::new(0.0, 1.5, 0.5)));
    }
}

#[test]
fn trimesh_plane_vertex_intersection() {
    let mesh = build_diamond(&Pose::identity());

    let result = mesh.intersection_with_local_plane(Vector::Z, 0.0, core::f32::EPSILON);

    assert!(matches!(result, IntersectResult::Intersect(_)));

    if let IntersectResult::Intersect(line) = result {
        // Need to check points individually since order is not guaranteed
        let vertices = line.vertices();
        assert_eq!(vertices.len(), 3);
        assert!(vertices.contains(&Vector::new(-2.0, -1.0, 0.0)));
        assert!(vertices.contains(&Vector::new(2.0, -1.0, 0.0)));
        assert!(vertices.contains(&Vector::new(0.0, 2.0, 0.0)));
    }
}

#[test]
fn trimesh_plane_mixed_intersection() {
    let mesh = build_diamond(&Pose::identity());

    let result = mesh.intersection_with_local_plane(Vector::X, 0.0, core::f32::EPSILON);

    assert!(matches!(result, IntersectResult::Intersect(_)));

    if let IntersectResult::Intersect(line) = result {
        // Need to check points individually since order is not guaranteed
        let vertices = line.vertices();
        assert_eq!(vertices.len(), 4);
        assert!(vertices.contains(&Vector::new(0.0, 2.0, 0.0)));
        assert!(vertices.contains(&Vector::new(0.0, 0.0, 2.0)));
        assert!(vertices.contains(&Vector::new(0.0, -1.0, 0.0)));
        assert!(vertices.contains(&Vector::new(0.0, 0.0, -2.0)));
    }
}

#[test]
fn trimesh_plane_multi_intersection() {
    let mut mesh = build_diamond(&Pose::identity());
    mesh.append(&build_diamond(&Pose::translation(-5.0, 0.0, 0.0)));

    let result = mesh.intersection_with_local_plane(Vector::Z, 0.5, core::f32::EPSILON);

    assert!(matches!(result, IntersectResult::Intersect(_)));

    if let IntersectResult::Intersect(line) = result {
        // Need to check points individually since order is not guaranteed
        let vertices = line.vertices();
        assert_eq!(vertices.len(), 6);

        assert!(vertices.contains(&Vector::new(-1.5, -0.75, 0.5)));
        assert!(vertices.contains(&Vector::new(1.5, -0.75, 0.5)));
        assert!(vertices.contains(&Vector::new(0.0, 1.5, 0.5)));

        assert!(vertices.contains(&Vector::new(-6.5, -0.75, 0.5)));
        assert!(vertices.contains(&Vector::new(-3.5, -0.75, 0.5)));
        assert!(vertices.contains(&Vector::new(-5.0, 1.5, 0.5)));
    }
}

#[test]
fn trimesh_plane_above() {
    let mesh = build_diamond(&Pose::identity());

    let result = mesh.intersection_with_local_plane(Vector::Z, -5.0, core::f32::EPSILON);

    assert!(matches!(result, IntersectResult::Positive));
}

#[test]
fn trimesh_plane_below() {
    let mesh = build_diamond(&Pose::identity());

    let result = mesh.intersection_with_local_plane(Vector::Z, 5.0, core::f32::EPSILON);

    assert!(matches!(result, IntersectResult::Negative));
}
