//! Point inclusion and projection.

pub use self::point_composite_shape::{
    PointCompositeShapeProjBestFirstVisitor, PointCompositeShapeProjWithFeatureBestFirstVisitor,
    PointCompositeShapeProjWithLocationBestFirstVisitor,
};
#[doc(inline)]
pub use self::point_query::{PointProjection, PointQuery, PointQueryWithLocation};
pub use self::point_support_map::local_point_projection_on_support_map;

mod point_aabb;
mod point_ball;
mod point_bounding_sphere;
mod point_capsule;
mod point_composite_shape;
#[cfg(feature = "dim3")]
mod point_cone;
mod point_cuboid;
#[cfg(feature = "dim3")]
mod point_cylinder;
mod point_halfspace;
mod point_heightfield;
#[doc(hidden)]
pub mod point_query;
mod point_round_shape;
mod point_segment;
mod point_support_map;
#[cfg(feature = "dim3")]
mod point_tetrahedron;
mod point_triangle;
