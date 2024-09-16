use crate::common::{generate, generate_trimesh_around_origin, unref};
use na::Isometry3;
use parry3d::bounding_volume::BoundingVolume;
use parry3d::bounding_volume::{Aabb, BoundingSphere};
use parry3d::shape::{
    Ball, Capsule, Cone, ConvexPolyhedron, Cuboid, Cylinder, Segment, TriMesh, Triangle,
};
use rand::SeedableRng;
use rand_isaac::IsaacRng;
use test::Bencher;

#[path = "../common/macros.rs"]
#[macro_use]
mod macros;

/*
 * Bounding volume methods.
 */
bench_method!(
    bench_aabb_intersects_aabb_always_true,
    intersects,
    aabb1: Aabb,
    aabb2: Aabb
);
bench_method!(
    bench_bounding_sphere_intersects_bounding_sphere_always_true,
    intersects,
    bs1: BoundingSphere,
    bs2: BoundingSphere
);

bench_method!(bench_aabb_contains_aabb, contains, aabb1: Aabb, aabb2: Aabb);
bench_method!(
    bench_bounding_sphere_contains_bounding_sphere,
    contains,
    bs1: BoundingSphere,
    bs2: BoundingSphere
);

bench_method!(bench_aabb_merged_aabb, merged, aabb1: Aabb, aabb2: Aabb);
bench_method!(
    bench_bounding_sphere_merged_bounding_sphere,
    merged,
    bs1: BoundingSphere,
    bs2: BoundingSphere
);

bench_method!(bench_aabb_loosened_aabb, loosened, aabb1: Aabb, margin: f32);
bench_method!(
    bench_bounding_sphere_loosened_bounding_sphere,
    loosened,
    bs1: BoundingSphere,
    margin: f32
);

/*
 * Bounding volume construction.
 */
bench_method!(bench_cuboid_aabb, aabb: Aabb, c: Cuboid, m: Isometry3<f32>);
bench_method!(
    bench_cuboid_bounding_sphere,
    bounding_sphere: BoundingSphere,
    c: Cuboid,
    m: Isometry3<f32>
);

bench_method!(bench_ball_aabb, aabb: Aabb, b: Ball, m: Isometry3<f32>);
bench_method!(
    bench_ball_bounding_sphere,
    bounding_sphere: BoundingSphere,
    b: Ball,
    m: Isometry3<f32>
);

bench_method!(
    bench_capsule_aabb,
    aabb: Aabb,
    c: Capsule,
    m: Isometry3<f32>
);
bench_method!(
    bench_capsule_bounding_sphere,
    bounding_sphere: BoundingSphere,
    c: Capsule,
    m: Isometry3<f32>
);

bench_method!(bench_cone_aabb, aabb: Aabb, c: Cone, m: Isometry3<f32>);
bench_method!(
    bench_cone_bounding_sphere,
    bounding_sphere: BoundingSphere,
    c: Cone,
    m: Isometry3<f32>
);

bench_method!(
    bench_cylinder_aabb,
    aabb: Aabb,
    c: Cylinder,
    m: Isometry3<f32>
);
bench_method!(
    bench_cylinder_bounding_sphere,
    bounding_sphere: BoundingSphere,
    c: Cylinder,
    m: Isometry3<f32>
);

bench_method!(
    bench_segment_aabb,
    aabb: Aabb,
    c: Segment,
    m: Isometry3<f32>
);
bench_method!(
    bench_segment_bounding_sphere,
    bounding_sphere: BoundingSphere,
    c: Segment,
    m: Isometry3<f32>
);

bench_method!(
    bench_triangle_aabb,
    aabb: Aabb,
    c: Triangle,
    m: Isometry3<f32>
);
bench_method!(
    bench_triangle_bounding_sphere,
    bounding_sphere: BoundingSphere,
    c: Triangle,
    m: Isometry3<f32>
);

bench_method!(
    bench_convex_aabb,
    aabb: Aabb,
    c: ConvexPolyhedron,
    m: Isometry3<f32>
);
bench_method!(
    bench_convex_bounding_sphere,
    bounding_sphere: BoundingSphere,
    c: ConvexPolyhedron,
    m: Isometry3<f32>
);

bench_method_gen!(
    bench_mesh_aabb,
    aabb: Aabb,
    mesh: TriMesh = generate_trimesh_around_origin,
    m: Isometry3<f32> = generate
);
bench_method_gen!(
    bench_mesh_bounding_sphere,
    bounding_sphere: BoundingSphere,
    mesh: TriMesh = generate_trimesh_around_origin,
    m: Isometry3<f32> = generate
);
