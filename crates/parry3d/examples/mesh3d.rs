use parry3d::math::Vector;
use parry3d::shape::TriMesh;

fn main() {
    let points = vec![
        Vector::new(0.0, 1.0, 0.0),
        Vector::new(-1.0, -0.5, 0.0),
        Vector::new(0.0, -0.5, -1.0),
        Vector::new(1.0, -0.5, 0.0),
    ];

    let indices = vec![[0u32, 1, 2], [0, 2, 3], [0, 3, 1]];

    // Build the mesh.
    let mesh = TriMesh::new(points, indices).unwrap();

    assert!(mesh.vertices().len() == 4);
}
