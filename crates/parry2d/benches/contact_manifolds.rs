// Criterion benchmarks for steady-state 2D contact-manifold updates,
// mirroring the parry3d suite: dedicated convex pairs, the generic PFM
// fallback, and trimesh-vs-convex persistent manifolds. The body is written
// against `parry2d::math::Real` so the parry2d-f64 crate reuses it verbatim.

use criterion::{criterion_group, criterion_main, Criterion};
use parry2d::math::{Pose, Real, Vector};
use parry2d::query::{
    ContactManifold, ContactManifoldsWorkspace, DefaultQueryDispatcher, PersistentQueryDispatcher,
};
use parry2d::shape::{Ball, Capsule, ConvexPolygon, Cuboid, TriMesh};
use std::hint::black_box;
use std::time::Duration;

type Manifold = ContactManifold<(), ()>;

fn drifting_pose(step: u32) -> Pose {
    let t = step as Real * 0.02;
    Pose::translation(t.sin() * 0.03, 1.2 + t.cos() * 0.02)
}

fn hexagon(radius: Real) -> ConvexPolygon {
    let pts: Vec<Vector> = (0..6)
        .map(|i| {
            let a = i as Real * core::f64::consts::PI as Real / 3.0;
            Vector::new(a.cos() * radius, a.sin() * radius)
        })
        .collect();
    ConvexPolygon::from_convex_polyline(pts).unwrap()
}

fn bench_convex_pairs(c: &mut Criterion) {
    let dispatcher = DefaultQueryDispatcher;
    let prediction = 0.1;

    let mut g = c.benchmark_group("contact_manifold_convex_convex");

    g.bench_function("cuboid_cuboid", |b| {
        let c1 = Cuboid::new(Vector::new(1.0, 1.0));
        let c2 = Cuboid::new(Vector::new(0.8, 0.8));
        let mut manifold = Manifold::new();
        let mut step = 0u32;
        b.iter(|| {
            step = step.wrapping_add(1);
            let pos12 = drifting_pose(step) * Pose::translation(0.0, 0.35);
            dispatcher
                .contact_manifold_convex_convex(
                    &pos12,
                    &c1,
                    &c2,
                    None,
                    None,
                    prediction,
                    &mut manifold,
                )
                .unwrap();
            black_box(manifold.points.len())
        })
    });

    // Hexagon/capsule: generic support-map fallback (GJK + EPA on penetration).
    g.bench_function("pfm_hexagon_capsule", |b| {
        let poly = hexagon(1.0);
        let cap = Capsule::new_y(0.6, 0.3);
        let mut manifold = Manifold::new();
        let mut step = 0u32;
        b.iter(|| {
            step = step.wrapping_add(1);
            let pos12 = drifting_pose(step);
            dispatcher
                .contact_manifold_convex_convex(
                    &pos12,
                    &poly,
                    &cap,
                    None,
                    None,
                    prediction,
                    &mut manifold,
                )
                .unwrap();
            black_box(manifold.points.len())
        })
    });

    g.bench_function("capsule_capsule", |b| {
        let cap1 = Capsule::new_y(0.8, 0.4);
        let cap2 = Capsule::new_y(0.6, 0.3);
        let mut manifold = Manifold::new();
        let mut step = 0u32;
        b.iter(|| {
            step = step.wrapping_add(1);
            let pos12 = drifting_pose(step) * Pose::translation(0.0, -0.2);
            dispatcher
                .contact_manifold_convex_convex(
                    &pos12,
                    &cap1,
                    &cap2,
                    None,
                    None,
                    prediction,
                    &mut manifold,
                )
                .unwrap();
            black_box(manifold.points.len())
        })
    });

    g.finish();
}

/// A triangulated band: a sine-curve top edge over a flat bottom edge.
fn make_terrain(columns: u32) -> TriMesh {
    let n = columns as usize;
    let mut vertices = Vec::with_capacity((n + 1) * 2);
    let mut indices = Vec::with_capacity(n * 2);
    for ix in 0..=n {
        let x = ix as Real / n as Real * 40.0 - 20.0;
        let y = (x * 0.8).sin() * 0.4;
        vertices.push(Vector::new(x, y)); // top: 2*ix
        vertices.push(Vector::new(x, -1.5)); // bottom: 2*ix + 1
    }
    for ix in 0..n {
        let top0 = (ix * 2) as u32;
        let bot0 = top0 + 1;
        let top1 = top0 + 2;
        let bot1 = top0 + 3;
        indices.push([top0, bot0, top1]);
        indices.push([bot0, bot1, top1]);
    }
    TriMesh::new(vertices, indices).unwrap()
}

fn bench_trimesh_shape(c: &mut Criterion) {
    let dispatcher = DefaultQueryDispatcher;
    let prediction = 0.1;
    let terrain = make_terrain(512); // 1024 triangles

    let mut g = c.benchmark_group("contact_manifolds_trimesh");

    g.bench_function("trimesh_vs_ball", |b| {
        let ball = Ball::new(0.7);
        let mut manifolds: Vec<Manifold> = Vec::new();
        let mut workspace: Option<ContactManifoldsWorkspace> = None;
        let mut step = 0u32;
        b.iter(|| {
            step = step.wrapping_add(1);
            let t = step as Real * 0.02;
            let pos12 = Pose::translation(t.sin() * 0.5, 0.55);
            dispatcher
                .contact_manifolds(
                    &pos12,
                    &terrain,
                    &ball,
                    prediction,
                    &mut manifolds,
                    &mut workspace,
                )
                .unwrap();
            black_box(manifolds.len())
        })
    });

    g.bench_function("trimesh_vs_cuboid", |b| {
        let cuboid = Cuboid::new(Vector::new(0.6, 0.6));
        let mut manifolds: Vec<Manifold> = Vec::new();
        let mut workspace: Option<ContactManifoldsWorkspace> = None;
        let mut step = 0u32;
        b.iter(|| {
            step = step.wrapping_add(1);
            let t = step as Real * 0.02;
            let pos12 = Pose::translation(t.sin() * 0.5, 0.5);
            dispatcher
                .contact_manifolds(
                    &pos12,
                    &terrain,
                    &cuboid,
                    prediction,
                    &mut manifolds,
                    &mut workspace,
                )
                .unwrap();
            black_box(manifolds.len())
        })
    });

    g.finish();
}

fn config() -> Criterion {
    Criterion::default()
        .warm_up_time(Duration::from_millis(800))
        .measurement_time(Duration::from_secs(3))
        .sample_size(60)
}

criterion_group! {
    name = benches;
    config = config();
    targets = bench_convex_pairs, bench_trimesh_shape
}
criterion_main!(benches);
