extern crate nalgebra as na;

use na::Point3;
use parry3d::shape::TriMesh;

fn main() {
    let points = vec![
        Point3::new(0.0, 1.0, 0.0),
        Point3::new(-1.0, -0.5, 0.0),
        Point3::new(0.0, -0.5, -1.0),
        Point3::new(1.0, -0.5, 0.0),
    ];

    let indices = vec![[0u32, 1, 2], [0, 2, 3], [0, 3, 1]];

    // Build the mesh.
    let mesh = TriMesh::new(points, indices).unwrap();

    assert!(mesh.vertices().len() == 4);
}
