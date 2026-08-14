// Criterion benchmarks for steady-state contact-manifold updates.
//
// These target the narrow-phase hot paths from the 2026-08 optimization audit:
// per-rebuild `manifold.points.clone()` (P1), per-contact `EPA::new()` (P2),
// and loop-invariant work in the trimesh manifold loop (P6). Poses drift a
// little every iteration so contact tracking stays on the realistic
// "rebuild + match old points" path.

use criterion::{criterion_group, criterion_main, Criterion};
use parry3d::math::{Pose, Real, Vector};
use parry3d::query::{
    ContactManifold, ContactManifoldsWorkspace, DefaultQueryDispatcher, PersistentQueryDispatcher,
};
use parry3d::shape::{Ball, Capsule, Cone, Cuboid, Cylinder, TriMesh};
use std::hint::black_box;
use std::time::Duration;

type Manifold = ContactManifold<(), ()>;

fn drifting_pose(step: u32) -> Pose {
    let t = step as Real * 0.02;
    Pose::translation(t.sin() * 0.03, 1.4 + t.cos() * 0.02, 0.0)
}

fn bench_convex_pairs(c: &mut Criterion) {
    let dispatcher = DefaultQueryDispatcher;
    let prediction = 0.1;

    let mut g = c.benchmark_group("contact_manifold_convex_convex");

    // Cuboid/cuboid: dedicated SAT-based path with the points.clone() pattern.
    g.bench_function("cuboid_cuboid", |b| {
        let c1 = Cuboid::new(Vector::new(1.0, 1.0, 1.0));
        let c2 = Cuboid::new(Vector::new(0.8, 0.8, 0.8));
        let mut manifold = Manifold::new();
        let mut step = 0u32;
        b.iter(|| {
            step = step.wrapping_add(1);
            let pos12 = drifting_pose(step) * Pose::translation(0.0, 0.35, 0.0);
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

    // Cylinder/cone: generic PFM fallback (GJK + EPA on penetration).
    g.bench_function("pfm_pfm_cylinder_cone", |b| {
        let c1 = Cylinder::new(0.8, 0.6);
        let c2 = Cone::new(0.7, 0.5);
        let mut manifold = Manifold::new();
        let mut step = 0u32;
        b.iter(|| {
            step = step.wrapping_add(1);
            let pos12 = drifting_pose(step);
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

    // Capsule/capsule: dedicated path, also carries the clone pattern.
    g.bench_function("capsule_capsule", |b| {
        let cap1 = Capsule::new_y(0.8, 0.4);
        let cap2 = Capsule::new_y(0.6, 0.3);
        let mut manifold = Manifold::new();
        let mut step = 0u32;
        b.iter(|| {
            step = step.wrapping_add(1);
            let pos12 = drifting_pose(step) * Pose::translation(0.0, -0.2, 0.0);
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

fn bench_trimesh_shape(c: &mut Criterion) {
    let dispatcher = DefaultQueryDispatcher;
    let prediction = 0.1;
    let terrain = make_terrain(48); // 4608 triangles

    // Same terrain with internal-edge fixing: pseudo-normals are computed and
    // fetched per candidate triangle every frame.
    let (fixed_vtx, fixed_idx) = {
        let t = make_terrain(48);
        (t.vertices().to_vec(), t.indices().to_vec())
    };
    let terrain_fixed = TriMesh::with_flags(
        fixed_vtx,
        fixed_idx,
        parry3d::shape::TriMeshFlags::FIX_INTERNAL_EDGES,
    )
    .unwrap();

    let mut g = c.benchmark_group("contact_manifolds_trimesh");

    g.bench_function("trimesh_fixed_edges_vs_cuboid", |b| {
        let cuboid = Cuboid::new(Vector::new(0.6, 0.6, 0.6));
        let mut manifolds: Vec<Manifold> = Vec::new();
        let mut workspace: Option<ContactManifoldsWorkspace> = None;
        let mut step = 0u32;
        b.iter(|| {
            step = step.wrapping_add(1);
            let t = step as Real * 0.02;
            let pos12 = Pose::translation(t.sin() * 0.5, 0.5, t.cos() * 0.5);
            dispatcher
                .contact_manifolds(
                    &pos12,
                    &terrain_fixed,
                    &cuboid,
                    prediction,
                    &mut manifolds,
                    &mut workspace,
                )
                .unwrap();
            black_box(manifolds.len())
        })
    });

    g.bench_function("trimesh_fixed_edges_vs_ball", |b| {
        let ball = Ball::new(0.7);
        let mut manifolds: Vec<Manifold> = Vec::new();
        let mut workspace: Option<ContactManifoldsWorkspace> = None;
        let mut step = 0u32;
        b.iter(|| {
            step = step.wrapping_add(1);
            let t = step as Real * 0.02;
            let pos12 = Pose::translation(t.sin() * 0.5, 0.55, t.cos() * 0.5);
            dispatcher
                .contact_manifolds(
                    &pos12,
                    &terrain_fixed,
                    &ball,
                    prediction,
                    &mut manifolds,
                    &mut workspace,
                )
                .unwrap();
            black_box(manifolds.len())
        })
    });

    g.bench_function("trimesh_vs_ball", |b| {
        let ball = Ball::new(0.7);
        let mut manifolds: Vec<Manifold> = Vec::new();
        let mut workspace: Option<ContactManifoldsWorkspace> = None;
        let mut step = 0u32;
        b.iter(|| {
            step = step.wrapping_add(1);
            let t = step as Real * 0.02;
            let pos12 = Pose::translation(t.sin() * 0.5, 0.55, t.cos() * 0.5);
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
        let cuboid = Cuboid::new(Vector::new(0.6, 0.6, 0.6));
        let mut manifolds: Vec<Manifold> = Vec::new();
        let mut workspace: Option<ContactManifoldsWorkspace> = None;
        let mut step = 0u32;
        b.iter(|| {
            step = step.wrapping_add(1);
            let t = step as Real * 0.02;
            let pos12 = Pose::translation(t.sin() * 0.5, 0.5, t.cos() * 0.5);
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
