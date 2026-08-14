// Criterion benchmarks for GJK/EPA contact queries between smooth convex
// shapes (support-map path). Deep-penetration cases exercise the per-call
// `EPA::new()` allocation (audit finding P2); separated/shallow cases stay
// in pure GJK (finding P5's per-iteration sqrt).

use criterion::{criterion_group, criterion_main, Criterion};
use parry3d::math::Pose;
use parry3d::query;
use parry3d::shape::{Capsule, Cone, Cylinder};
use std::hint::black_box;
use std::time::Duration;

fn bench_gjk_epa(c: &mut Criterion) {
    let cyl = Cylinder::new(0.8, 0.6);
    let cone = Cone::new(0.7, 0.5);
    let cap = Capsule::new_y(0.8, 0.4);
    let identity = Pose::identity();

    let mut g = c.benchmark_group("gjk_epa_contact");

    g.bench_function("gjk_separated", |b| {
        let pos2 = Pose::translation(2.6, 0.4, 0.1);
        b.iter(|| black_box(query::contact(&identity, &cyl, &pos2, &cone, 1.5).unwrap()))
    });

    g.bench_function("gjk_shallow_penetration", |b| {
        let pos2 = Pose::translation(1.02, 0.05, 0.0);
        b.iter(|| black_box(query::contact(&identity, &cyl, &pos2, &cone, 0.1).unwrap()))
    });

    g.bench_function("epa_deep_penetration", |b| {
        let pos2 = Pose::translation(0.35, 0.1, 0.05);
        b.iter(|| black_box(query::contact(&identity, &cyl, &pos2, &cone, 0.1).unwrap()))
    });

    g.bench_function("epa_deep_capsule_cylinder", |b| {
        let pos2 = Pose::translation(0.25, 0.2, 0.0);
        b.iter(|| black_box(query::contact(&identity, &cap, &pos2, &cyl, 0.1).unwrap()))
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
