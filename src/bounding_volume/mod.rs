//! Bounding volumes.

#[doc(inline)]
pub use crate::bounding_volume::aabb::AABB;
pub use crate::bounding_volume::simd_aabb::SimdAABB;

#[doc(inline)]
pub use crate::bounding_volume::bounding_sphere::BoundingSphere;
#[doc(inline)]
pub use crate::bounding_volume::bounding_volume::BoundingVolume;

#[doc(hidden)]
pub mod bounding_volume;

#[doc(hidden)]
pub mod aabb;
mod aabb_ball;
#[cfg(feature = "dim2")]
mod aabb_convex_polygon;
#[cfg(feature = "dim3")]
mod aabb_convex_polyhedron;
mod aabb_cuboid;
mod aabb_halfspace;
mod aabb_heightfield;
mod aabb_support_map;
mod aabb_triangle;
mod aabb_utils;

mod aabb_capsule;
#[doc(hidden)]
pub mod bounding_sphere;
mod bounding_sphere_ball;
mod bounding_sphere_capsule;
#[cfg(feature = "dim3")]
mod bounding_sphere_cone;
#[cfg(feature = "dim3")]
mod bounding_sphere_convex;
#[cfg(feature = "dim2")]
mod bounding_sphere_convex_polygon;
mod bounding_sphere_cuboid;
#[cfg(feature = "dim3")]
mod bounding_sphere_cylinder;
mod bounding_sphere_halfspace;
mod bounding_sphere_heightfield;
mod bounding_sphere_polyline;
mod bounding_sphere_segment;
mod bounding_sphere_triangle;
mod bounding_sphere_trimesh;
mod bounding_sphere_utils;
mod simd_aabb;

/// Free functions for some special cases of bounding-volume computation.
pub mod details {
    #[cfg(feature = "dim3")]
    pub use super::aabb_utils::support_map_aabb;
    pub use super::aabb_utils::{local_point_cloud_aabb, local_support_map_aabb, point_cloud_aabb};
    pub use super::bounding_sphere_utils::point_cloud_bounding_sphere;
}
