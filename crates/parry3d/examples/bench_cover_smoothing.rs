use parry3d::shape::Ball;
use parry3d::transformation::{volume_mesh, MeshEnclosure, VolumeMeshParameters};
use std::time::Instant;

fn main() {
    let (vertices, indices) = Ball::new(1.0).to_trimesh(30, 30);
    let mut params = VolumeMeshParameters::new(0.2);
    params.enclosure = MeshEnclosure::Cover;
    params.cover_subdivisions = 2;

    let start = Instant::now();
    let raw = volume_mesh(&vertices, &indices, &params).unwrap();
    let _warmup = volume_mesh(&vertices, &indices, &params).unwrap();
    let start = Instant::now();
    let raw2 = volume_mesh(&vertices, &indices, &params).unwrap();
    let build = start.elapsed();
    let _ = raw2;

    params.cover_smoothing = 30;
    let start = Instant::now();
    let smoothed = volume_mesh(&vertices, &indices, &params).unwrap();
    let total = start.elapsed();

    let checksum: f64 = smoothed
        .vertices
        .iter()
        .map(|v| (v.x as f64) + (v.y as f64) * 3.0 + (v.z as f64) * 7.0)
        .sum();
    println!(
        "{} cells; build {build:?}, build+smooth {total:?}, smoothing {:?}, checksum {checksum:.9}",
        smoothed.cells.len(),
        total - build
    );
    let _ = raw;
}
