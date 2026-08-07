// Regression test for https://github.com/dimforge/parry/issues/177
//
// `Debug` for `TriMesh` used to print just "GenericTriMesh"; it now summarizes the
// vertex/triangle counts, local AABB, flags, and optional topology data.

use parry3d::shape::{TriMesh, TriMeshFlags};

use parry3d::math::Vector;

#[test]
fn trimesh_debug_prints_summary() {
    let vertices = vec![
        Vector::new(0.0, 0.0, 0.0),
        Vector::new(1.0, 0.0, 0.0),
        Vector::new(0.0, 1.0, 0.0),
        Vector::new(1.0, 1.0, 0.0),
    ];
    let indices = vec![[0, 1, 2], [1, 3, 2]];
    let mesh = TriMesh::with_flags(vertices, indices, TriMeshFlags::HALF_EDGE_TOPOLOGY).unwrap();

    let dbg = format!("{mesh:?}");
    assert!(dbg.contains("TriMesh"), "{dbg}");
    assert!(dbg.contains("num_vertices: 4"), "{dbg}");
    assert!(dbg.contains("num_triangles: 2"), "{dbg}");
    assert!(dbg.contains("local_aabb"), "{dbg}");
    assert!(dbg.contains("flags"), "{dbg}");
    assert!(dbg.contains("has_pseudo_normals: false"), "{dbg}");
    assert!(dbg.contains("has_topology: true"), "{dbg}");
    assert!(dbg.contains("has_connected_components: false"), "{dbg}");
}
