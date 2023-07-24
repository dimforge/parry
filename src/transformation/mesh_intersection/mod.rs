pub use self::mesh_intersection::intersect_meshes;
pub use self::mesh_intersection_error::MeshIntersectionError;
pub(self) use triangle_triangle_intersection::*;

use crate::math::{Real, real};

mod mesh_intersection;
mod mesh_intersection_error;
mod triangle_triangle_intersection;

pub(self) const EPS: Real = real!(1.0e-6);
