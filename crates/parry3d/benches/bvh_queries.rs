// Criterion benchmarks for BVH-backed queries on a trimesh: AABB
// interference enumeration, ray casts, and point projection. These exercise
// the BvhNode predicate hot paths (intersects/contains/cast_ray) through the
// public TriMesh query API. Written against `Real` so f64 crates can reuse
// the body through an alias shim.

use criterion::{criterion_group, criterion_main, Criterion};
use parry3d::bounding_volume::Aabb;
use parry3d::math::{Real, Vector};
use parry3d::query::{PointQuery, Ray, RayCast};
use parry3d::shape::TriMesh;
use std::hint::black_box;
use std::time::Duration;

fn make_terrain(subdivisions: u32) -> TriMesh {
    let n = subdivisions as usize;
    let mut vertices = Vec::with_capacity((n + 1) * (n + 1));
    let mut indices = Vec::with_capacity(n * n * 2);
    for iz in 0..=n {
        for ix in 0..=n {
            let x = ix as Real / n as Real * 20.0 - 10.0;
            let z = iz as Real / n as Real * 20.0 - 10.0;
            let y = (x * 0.8).sin() * 0.4 + (z * 0.7).cos() * 0.4;
            vertices.push(Vector::new(x, y, z));
        }
    }
    let stride = n + 1;
    for iz in 0..n {
        for ix in 0..n {
            let a = (iz * stride + ix) as u32;
            let b = (iz * stride + ix + 1) as u32;
            let c = ((iz + 1) * stride + ix) as u32;
            let d = ((iz + 1) * stride + ix + 1) as u32;
            indices.push([a, b, c]);
            indices.push([b, d, c]);
        }
    }
    TriMesh::new(vertices, indices).unwrap()
}

fn bench_bvh_queries(c: &mut Criterion) {
    let terrain = make_terrain(64); // 8192 triangles
    let bvh = terrain.bvh();

    let mut g = c.benchmark_group("bvh_queries");

    g.bench_function("intersect_aabb_64", |b| {
        let queries: Vec<Aabb> = (0..64)
            .map(|i| {
                let t = i as Real * 0.37;
                let center = Vector::new(t.sin() * 8.0, 0.0, t.cos() * 8.0);
                let he = Vector::new(0.8, 2.0, 0.8);
                Aabb::new(center - he, center + he)
            })
            .collect();
        b.iter(|| {
            let mut count = 0usize;
            for q in &queries {
                count += bvh.intersect_aabb(q).count();
            }
            black_box(count)
        })
    });

    g.bench_function("cast_ray_64", |b| {
        let rays: Vec<Ray> = (0..64)
            .map(|i| {
                let t = i as Real * 0.61;
                Ray::new(
                    Vector::new(t.sin() * 9.0, 5.0, t.cos() * 9.0),
                    Vector::new(t.cos() * 0.1, -1.0, t.sin() * 0.1),
                )
            })
            .collect();
        b.iter(|| {
            let mut acc = 0.0;
            for r in &rays {
                if let Some(toi) = terrain.cast_local_ray(r, 100.0, true) {
                    acc += toi;
                }
            }
            black_box(acc)
        })
    });

    g.bench_function("project_point_64", |b| {
        let points: Vec<Vector> = (0..64)
            .map(|i| {
                let t = i as Real * 0.53;
                Vector::new(t.sin() * 8.0, 1.5 + t.cos(), t.cos() * 7.0)
            })
            .collect();
        b.iter(|| {
            let mut acc = 0.0;
            for p in &points {
                let proj = terrain.project_local_point(*p, true);
                acc += proj.point.x;
            }
            black_box(acc)
        })
    });

    g.finish();
}

fn config() -> Criterion {
    Criterion::default()
        .warm_up_time(Duration::from_millis(600))
        .measurement_time(Duration::from_secs(3))
        .sample_size(60)
}

criterion_group! {
    name = benches;
    config = config();
    targets = bench_bvh_queries
}
criterion_main!(benches);
