// Criterion benchmarks for 2D GJK/EPA contact queries between smooth convex
// shapes (support-map path), mirroring the parry3d suite. Written against
// `parry2d::math::Real` so parry2d-f64 reuses the body verbatim.

use criterion::{criterion_group, criterion_main, Criterion};
use parry2d::math::{Pose, Real, Vector};
use parry2d::query;
use parry2d::shape::{Capsule, ConvexPolygon};
use std::hint::black_box;
use std::time::Duration;

fn hexagon(radius: Real) -> ConvexPolygon {
    let pts: Vec<Vector> = (0..6)
        .map(|i| {
            let a = i as Real * core::f64::consts::PI as Real / 3.0;
            Vector::new(a.cos() * radius, a.sin() * radius)
        })
        .collect();
    ConvexPolygon::from_convex_polyline(pts).unwrap()
}

fn bench_gjk_epa(c: &mut Criterion) {
    let poly = hexagon(1.0);
    let cap = Capsule::new_y(0.8, 0.4);
    let identity = Pose::identity();

    let mut g = c.benchmark_group("gjk_epa_contact");

    g.bench_function("gjk_separated", |b| {
        let pos2 = Pose::translation(2.9, 0.4);
        b.iter(|| black_box(query::contact(&identity, &poly, &pos2, &cap, 1.5).unwrap()))
    });

    g.bench_function("gjk_shallow_penetration", |b| {
        let pos2 = Pose::translation(1.32, 0.05);
        b.iter(|| black_box(query::contact(&identity, &poly, &pos2, &cap, 0.1).unwrap()))
    });

    g.bench_function("epa_deep_penetration", |b| {
        let pos2 = Pose::translation(0.35, 0.1);
        b.iter(|| black_box(query::contact(&identity, &poly, &pos2, &cap, 0.1).unwrap()))
    });

    g.finish();
}

fn config() -> Criterion {
    Criterion::default()
        .warm_up_time(Duration::from_millis(500))
        .measurement_time(Duration::from_secs(2))
        .sample_size(80)
}

criterion_group! {
    name = benches;
    config = config();
    targets = bench_gjk_epa
}
criterion_main!(benches);
