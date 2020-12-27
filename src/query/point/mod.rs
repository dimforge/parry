//! Point inclusion and projection.

#[doc(inline)]
pub use self::point_query::{PointProjection, PointQuery, PointQueryWithLocation};
pub use self::point_support_map::local_point_projection_on_support_map;

mod point_aabb;
mod point_ball;
mod point_bounding_sphere;
mod point_capsule;
mod point_cuboid;
mod point_halfspace;
mod point_heightfield;
mod point_polyline;
#[doc(hidden)]
pub mod point_query;
mod point_round_shape;
mod point_segment;
mod point_support_map;
#[cfg(feature = "dim3")]
mod point_tetrahedron;
mod point_triangle;
mod point_trimesh;
