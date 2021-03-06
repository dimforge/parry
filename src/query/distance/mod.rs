//! Implementation details of the `distance` function.

pub use self::distance::distance;
pub use self::distance_ball_ball::distance_ball_ball;
pub use self::distance_composite_shape_shape::{
    distance_composite_shape_shape, distance_shape_composite_shape,
    CompositeShapeAgainstAnyDistanceVisitor,
};
pub use self::distance_halfspace_support_map::{
    distance_halfspace_support_map, distance_support_map_halfspace,
};
pub use self::distance_support_map_support_map::{
    distance_support_map_support_map, distance_support_map_support_map_with_params,
};

mod distance;
mod distance_ball_ball;
mod distance_composite_shape_shape;
mod distance_halfspace_support_map;
mod distance_support_map_support_map;
