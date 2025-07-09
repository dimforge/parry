use parry3d::math::Point;
use parry3d::shape::{TriMesh, TriMeshFlags};

#[test]
// From https://github.com/dimforge/parry/issues/115
fn mesh_connected_components_grouped_faces() {
    let verts = vec![
        // Face 0
        Point::new(15.82, 6.455, -0.15), // <- Vertex shared with face 1.
        Point::new(9.915, 6.455, -0.15),
        Point::new(9.915, 6.4, 0.0), // <- Vertex shared with face 1.
        // Face1
        Point::new(15.82, 6.455, -0.15), // <- Vertex shared with face 0.
        Point::new(9.915, 6.4, 0.0),     // <- Vertex shared with face 0.
        Point::new(15.82, 6.4, 0.0),
    ];

    let mut roof = TriMesh::new(verts, vec![[0, 1, 2], [3, 4, 5]]).unwrap();

    if let Err(e) =
        roof.set_flags(TriMeshFlags::MERGE_DUPLICATE_VERTICES | TriMeshFlags::CONNECTED_COMPONENTS)
    {
        panic!("{:?}", e);
    }

    let components = roof.connected_components().unwrap();
    println!("components: {components:?}");
    assert_eq!(components.ranges.len(), 2); // Only one connected-component (two range values).
    assert_eq!(components.grouped_faces.len(), 2); // Only two faces in the connected-component.
}
