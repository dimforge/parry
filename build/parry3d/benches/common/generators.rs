use na::{Point2, Point3};
use parry3d::shape::TriMesh;
use rand::Rng;
use rand_isaac::IsaacRng;

pub fn generate_trimesh_around_origin<R: Rng>(rng: &mut R) -> TriMesh {
    let count = 1000;
    let pts = (0..3 * count).map(|_| rng.gen::<Point3<f32>>() * 3.0).collect();
    let indices = (0..count)
        .map(|i| [i * 3, i * 3 + 1, i * 3 + 2])
        .collect();

    TriMesh::new(pts, indices)
}
pub fn generate_trimesh_around_origin_100<R: Rng>(rng: &mut R) -> TriMesh {
    let count = 100;
    let pts = (0..3 * count).map(|_| rng.gen::<Point3<f32>>() * 3.0).collect();
    let indices = (0..count)
        .map(|i| [i * 3, i * 3 + 1, i * 3 + 2])
        .collect();

    TriMesh::new(pts, indices)
}

pub fn generate_trimesh_around_origin_10000<R: Rng>(rng: &mut R) -> TriMesh {
    let count = 10000;
    let pts = (0..3 * count).map(|_| rng.gen::<Point3<f32>>() * 3.0).collect();
    let indices = (0..count)
        .map(|i| [i * 3, i * 3 + 1, i * 3 + 2])
        .collect();

    TriMesh::new(pts, indices)
}
