use na::{Point3, Vector3};
use parry3d::math::{Isometry, Real};
use parry3d::query::IntersectResult;
use parry3d::shape::TriMesh;

fn build_diamond(position: &Isometry<Real>) -> TriMesh {
    // Two tetrahedrons sharing a face
    let points = vec![
        position * Point3::new(0.0, 2.0, 0.0),
        position * Point3::new(-2.0, -1.0, 0.0),
        position * Point3::new(0.0, 0.0, 2.0),
        position * Point3::new(2.0, -1.0, 0.0),
        position * Point3::new(0.0, 0.0, -2.0),
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
    let mesh = build_diamond(&Isometry::identity());

    let result = mesh.intersection_with_local_plane(&Vector3::ith_axis(2), 0.5, std::f32::EPSILON);

    assert!(matches!(result, IntersectResult::Intersect(_)));

    if let IntersectResult::Intersect(line) = result {
        // Need to check points individually since order is not guaranteed
        let vertices = line.vertices();
        assert_eq!(vertices.len(), 3);
        assert!(vertices.contains(&Point3::new(-1.5, -0.75, 0.5)));
        assert!(vertices.contains(&Point3::new(1.5, -0.75, 0.5)));
        assert!(vertices.contains(&Point3::new(0.0, 1.5, 0.5)));
    }
}

#[test]
fn trimesh_plane_vertex_intersection() {
    let mesh = build_diamond(&Isometry::identity());

    let result = mesh.intersection_with_local_plane(&Vector3::ith_axis(2), 0.0, std::f32::EPSILON);

    assert!(matches!(result, IntersectResult::Intersect(_)));

    if let IntersectResult::Intersect(line) = result {
        // Need to check points individually since order is not guaranteed
        let vertices = line.vertices();
        assert_eq!(vertices.len(), 3);
        assert!(vertices.contains(&Point3::new(-2.0, -1.0, 0.0)));
        assert!(vertices.contains(&Point3::new(2.0, -1.0, 0.0)));
        assert!(vertices.contains(&Point3::new(0.0, 2.0, 0.0)));
    }
}

#[test]
fn trimesh_plane_mixed_intersection() {
    let mesh = build_diamond(&Isometry::identity());

    let result = mesh.intersection_with_local_plane(&Vector3::ith_axis(0), 0.0, std::f32::EPSILON);

    assert!(matches!(result, IntersectResult::Intersect(_)));

    if let IntersectResult::Intersect(line) = result {
        // Need to check points individually since order is not guaranteed
        let vertices = line.vertices();
        assert_eq!(vertices.len(), 4);
        assert!(vertices.contains(&Point3::new(0.0, 2.0, 0.0)));
        assert!(vertices.contains(&Point3::new(0.0, 0.0, 2.0)));
        assert!(vertices.contains(&Point3::new(0.0, -1.0, 0.0)));
        assert!(vertices.contains(&Point3::new(0.0, 0.0, -2.0)));
    }
}

#[test]
fn trimesh_plane_multi_intersection() {
    let mut mesh = build_diamond(&Isometry::identity());
    mesh.append(&build_diamond(&Isometry::translation(-5.0, 0.0, 0.0)));

    let result = mesh.intersection_with_local_plane(&Vector3::ith_axis(2), 0.5, std::f32::EPSILON);

    assert!(matches!(result, IntersectResult::Intersect(_)));

    if let IntersectResult::Intersect(line) = result {
        // Need to check points individually since order is not guaranteed
        let vertices = line.vertices();
        assert_eq!(vertices.len(), 6);

        assert!(vertices.contains(&Point3::new(-1.5, -0.75, 0.5)));
        assert!(vertices.contains(&Point3::new(1.5, -0.75, 0.5)));
        assert!(vertices.contains(&Point3::new(0.0, 1.5, 0.5)));

        assert!(vertices.contains(&Point3::new(-6.5, -0.75, 0.5)));
        assert!(vertices.contains(&Point3::new(-3.5, -0.75, 0.5)));
        assert!(vertices.contains(&Point3::new(-5.0, 1.5, 0.5)));
    }
}

#[test]
fn trimesh_plane_above() {
    let mesh = build_diamond(&Isometry::identity());

    let result = mesh.intersection_with_local_plane(&Vector3::ith_axis(2), -5.0, std::f32::EPSILON);

    assert!(matches!(result, IntersectResult::Positive));
}

#[test]
fn trimesh_plane_below() {
    let mesh = build_diamond(&Isometry::identity());

    let result = mesh.intersection_with_local_plane(&Vector3::ith_axis(2), 5.0, std::f32::EPSILON);

    assert!(matches!(result, IntersectResult::Negative));
}
