use crate::common::{generate, unref};
use parry3d::math::{Pose3, Vec3};
use parry3d::shape::ConvexPolyhedron;
use parry3d::shape::SupportMap;
use parry3d::shape::{Ball, Capsule, Cone, Cuboid, Cylinder, Segment, Triangle};
use rand::SeedableRng;
use rand_isaac::IsaacRng;
use test::Bencher;

#[path = "../common/macros.rs"]
#[macro_use]
mod macros;

bench_method!(
    bench_ball_support_map,
    support_point,
    c: Ball,
    m: Pose3,
    dir: Vec3
);
bench_method!(
    bench_cuboid_support_map,
    support_point,
    c: Cuboid,
    m: Pose3,
    dir: Vec3
);
bench_method!(
    bench_capsule_support_map,
    support_point,
    c: Capsule,
    m: Pose3,
    dir: Vec3
);
bench_method!(
    bench_cone_support_map,
    support_point,
    c: Cone,
    m: Pose3,
    dir: Vec3
);
bench_method!(
    bench_cylinder_support_map,
    support_point,
    c: Cylinder,
    m: Pose3,
    dir: Vec3
);
bench_method!(
    bench_segment_support_map,
    support_point,
    c: Segment,
    m: Pose3,
    dir: Vec3
);
bench_method!(
    bench_triangle_support_map,
    support_point,
    c: Triangle,
    m: Pose3,
    dir: Vec3
);
bench_method!(
    bench_convex_support_map,
    support_point,
    c: ConvexPolyhedron,
    m: Pose3,
    dir: Vec3
);
