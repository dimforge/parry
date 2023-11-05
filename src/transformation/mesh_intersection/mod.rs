use crate::math::Real;

// std required mostly due to spade which is not no_std

#[cfg(any(feature = "std", feature = "alloc"))]
pub use self::mesh_intersection::intersect_meshes;
pub use self::mesh_intersection_error::MeshIntersectionError;
pub(self) use triangle_triangle_intersection::*;

#[cfg(any(feature = "std", feature = "alloc"))]
mod mesh_intersection;
mod mesh_intersection_error;

mod triangle_triangle_intersection;

pub(self) const EPS: Real = 1.0e-6;
