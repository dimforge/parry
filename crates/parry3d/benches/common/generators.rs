use na::Point3;
use parry3d::shape::TriMesh;
use rand::Rng;

pub fn generate_trimesh_around_origin<R: Rng>(rng: &mut R) -> TriMesh {
    let pts = (0..3000)
        .map(|_| rng.random::<Point3<f32>>() * 3.0)
        .collect();
    let indices = (0..1000).map(|i| [i * 3, i * 3 + 1, i * 3 + 2]).collect();

    TriMesh::new(pts, indices).unwrap()
}
