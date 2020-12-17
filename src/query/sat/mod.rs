//! Application of the Separating-Axis-Theorem (SAT).

pub use self::sat_cuboid_cuboid::*;
pub use self::sat_cuboid_point::*;
pub use self::sat_cuboid_segment::*;
pub use self::sat_cuboid_support_map::*;
pub use self::sat_cuboid_triangle::*;
pub use self::sat_support_map_support_map::*;
#[cfg(feature = "dim3")]
pub use self::sat_triangle_segment::*;
// pub use self::sat_polygon_polygon::*;

mod sat_cuboid_cuboid;
mod sat_cuboid_point;
mod sat_cuboid_segment;
mod sat_cuboid_support_map;
mod sat_cuboid_triangle;
mod sat_support_map_support_map;
#[cfg(feature = "dim3")]
mod sat_triangle_segment;
// mod sat_polygon_polygon;
