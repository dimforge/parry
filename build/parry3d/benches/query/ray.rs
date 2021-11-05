use crate::common::{generate, generate_trimesh_around_origin,generate_trimesh_around_origin_100, generate_trimesh_around_origin_10000, unref};
use na::Isometry3;
use parry3d::bounding_volume::{BoundingSphere, AABB};
use parry3d::query::{Ray, RayCast};
use parry3d::shape::{
    Ball, Capsule, Cone, ConvexPolyhedron, Cuboid, Cylinder, Segment, TriMesh, Triangle,
};
use rand::SeedableRng;
use rand_isaac::IsaacRng;
use test::Bencher;

#[path = "../common/macros.rs"]
#[macro_use]
mod macros;

// FIXME: will the randomness of `solid` and `max_toi` affect too much the benchmark?
bench_method!(
    bench_ray_against_ball,
    cast_ray,
    b: Ball,
    pos: Isometry3<f32>,
    ray: Ray,
    max_toi: f32,
    solid: bool
);

bench_method!(
    bench_ray_against_cuboid,
    cast_ray,
    c: Cuboid,
    pos: Isometry3<f32>,
    ray: Ray,
    max_toi: f32,
    solid: bool
);

bench_method!(
    bench_ray_against_capsule,
    cast_ray,
    c: Capsule,
    pos: Isometry3<f32>,
    ray: Ray,
    max_toi: f32,
    solid: bool
);

bench_method!(
    bench_ray_against_cone,
    cast_ray,
    c: Cone,
    pos: Isometry3<f32>,
    ray: Ray,
    max_toi: f32,
    solid: bool
);

bench_method!(
    bench_ray_against_cylinder,
    cast_ray,
    c: Cylinder,
    pos: Isometry3<f32>,
    ray: Ray,
    max_toi: f32,
    solid: bool
);

bench_method!(
    bench_ray_against_aabb,
    cast_ray,
    a: AABB,
    pos: Isometry3<f32>,
    ray: Ray,
    max_toi: f32,
    solid: bool
);

bench_method!(
    bench_ray_against_bounding_sphere,
    cast_ray,
    b: BoundingSphere,
    pos: Isometry3<f32>,
    ray: Ray,
    max_toi: f32,
    solid: bool
);

bench_method!(
    bench_ray_against_ball_with_normal,
    cast_ray_and_get_normal,
    b: Ball,
    pos: Isometry3<f32>,
    ray: Ray,
    max_toi: f32,
    solid: bool
);

bench_method!(
    bench_ray_against_cuboid_with_normal,
    cast_ray_and_get_normal,
    c: Cuboid,
    pos: Isometry3<f32>,
    ray: Ray,
    max_toi: f32,
    solid: bool
);

bench_method!(
    bench_ray_against_capsule_with_normal,
    cast_ray_and_get_normal,
    c: Capsule,
    pos: Isometry3<f32>,
    ray: Ray,
    max_toi: f32,
    solid: bool
);

bench_method!(
    bench_ray_against_cone_with_normal,
    cast_ray_and_get_normal,
    c: Cone,
    pos: Isometry3<f32>,
    ray: Ray,
    max_toi: f32,
    solid: bool
);

bench_method!(
    bench_ray_against_cylinder_with_normal,
    cast_ray_and_get_normal,
    c: Cylinder,
    pos: Isometry3<f32>,
    ray: Ray,
    max_toi: f32,
    solid: bool
);

bench_method!(
    bench_ray_against_segment_with_normal,
    cast_ray_and_get_normal,
    c: Segment,
    pos: Isometry3<f32>,
    ray: Ray,
    max_toi: f32,
    solid: bool
);

bench_method!(
    bench_ray_against_triangle_with_normal,
    cast_ray_and_get_normal,
    c: Triangle,
    pos: Isometry3<f32>,
    ray: Ray,
    max_toi: f32,
    solid: bool
);

bench_method!(
    bench_ray_against_convex_with_normal,
    cast_ray_and_get_normal,
    c: ConvexPolyhedron,
    pos: Isometry3<f32>,
    ray: Ray,
    max_toi: f32,
    solid: bool
);

bench_method_gen!(
    bench_ray_against_trimesh_with_normal,
    cast_ray_and_get_normal,
    m: TriMesh = generate_trimesh_around_origin,
    pos: Isometry3<f32> = generate,
    ray: Ray = generate,
    max_toi: f32 = generate,
    solid: bool = generate
);

bench_method_gen!(
    bench_ray_against_trimesh_with_normal_100,
    cast_ray_and_get_normal,
    m: TriMesh = generate_trimesh_around_origin_100,
    pos: Isometry3<f32> = generate,
    ray: Ray = generate,
    max_toi: f32 = generate,
    solid: bool = generate
);

bench_method_gen!(
    bench_ray_against_trimesh_with_normal_10000,
    cast_ray_and_get_normal,
    m: TriMesh = generate_trimesh_around_origin_10000,
    pos: Isometry3<f32> = generate,
    ray: Ray = generate,
    max_toi: f32 = generate,
    solid: bool = generate
);
