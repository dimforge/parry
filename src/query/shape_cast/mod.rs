//! Implementation details of the `cast_shapes` function.

pub use self::shape_cast::{
    cast_shapes, cast_shapes_with_dispatcher, ShapeCastHit, ShapeCastOptions, ShapeCastStatus,
};
pub use self::shape_cast_ball_ball::cast_shapes_ball_ball;
pub use self::shape_cast_halfspace_support_map::{
    cast_shapes_halfspace_support_map, cast_shapes_support_map_halfspace,
};
#[cfg(feature = "std")]
pub use self::{
    shape_cast_composite_shape_shape::{
        cast_shapes_composite_shape_shape, cast_shapes_shape_composite_shape,
        TOICompositeShapeShapeBestFirstVisitor,
    },
    shape_cast_heightfield_shape::{cast_shapes_heightfield_shape, cast_shapes_shape_heightfield},
    shape_cast_support_map_support_map::cast_shapes_support_map_support_map,
};

mod shape_cast;
mod shape_cast_ball_ball;
#[cfg(feature = "std")]
mod shape_cast_composite_shape_shape;
mod shape_cast_halfspace_support_map;
#[cfg(feature = "std")]
mod shape_cast_heightfield_shape;
#[cfg(feature = "std")]
mod shape_cast_support_map_support_map;
