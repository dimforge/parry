pub use self::mesh_intersection::{intersect_meshes, tracked_intersection};
pub use self::mesh_intersection_error::MeshIntersectionError;
pub(self) use triangle_triangle_intersection::*;

use crate::math::Real;

mod mesh_intersection;
mod mesh_intersection_error;
mod triangle_triangle_intersection;

pub(self) const EPS: Real = 1.0e-6;
