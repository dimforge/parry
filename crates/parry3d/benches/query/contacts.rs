use crate::common::{generate, unref};
use na::Isometry3;
use parry3d::query;
use parry3d::shape::{Ball, Capsule, Cone, Cuboid, Cylinder};
use rand::SeedableRng;
use rand_isaac::IsaacRng;
use test::Bencher;

#[path = "../common/macros.rs"]
#[macro_use]
mod macros;

bench_free_fn!(
    bench_ball_against_ball,
    query::contact,
    pos1: Isometry3<f32>,
    b1: Ball,
    pos2: Isometry3<f32>,
    b2: Ball,
    prediction: f32
);

bench_free_fn!(
    bench_cuboid_against_cuboid,
    query::contact,
    pos1: Isometry3<f32>,
    b1: Cuboid,
    pos2: Isometry3<f32>,
    b2: Cuboid,
    prediction: f32
);

bench_free_fn!(
    bench_capsule_against_capsule,
    query::contact,
    pos1: Isometry3<f32>,
    b1: Capsule,
    pos2: Isometry3<f32>,
    b2: Capsule,
    prediction: f32
);

bench_free_fn!(
    bench_cone_against_cone,
    query::contact,
    pos1: Isometry3<f32>,
    b1: Cone,
    pos2: Isometry3<f32>,
    b2: Cone,
    prediction: f32
);

bench_free_fn!(
    bench_cylinder_against_cylinder,
    query::contact,
    pos1: Isometry3<f32>,
    b1: Cylinder,
    pos2: Isometry3<f32>,
    b2: Cylinder,
    prediction: f32
);
