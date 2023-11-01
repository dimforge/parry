use crate::math::Real;

// std required mostly due to spade which is not no_std

#[cfg(feature = "std")]
pub use self::mesh_intersection::intersect_meshes;
#[cfg(feature = "std")]
pub use self::mesh_intersection_error::MeshIntersectionError;
pub(self) use triangle_triangle_intersection::*;

#[cfg(feature = "std")]
mod mesh_intersection;
#[cfg(feature = "std")]
mod mesh_intersection_error;

mod triangle_triangle_intersection;

pub(self) const EPS: Real = 1.0e-6;
