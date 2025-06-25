//! Implementation details of the `cast_shapes_nonlinear` function.

#[cfg(feature = "alloc")]
pub use self::nonlinear_shape_cast_composite_shape_shape::{
    cast_shapes_nonlinear_composite_shape_shape, cast_shapes_nonlinear_shape_composite_shape,
    NonlinearTOICompositeShapeShapeBestFirstVisitor,
};
#[cfg(feature = "alloc")]
pub use self::nonlinear_shape_cast_voxels_shape::{
    cast_shapes_nonlinear_shape_voxels, cast_shapes_nonlinear_voxels_shape,
};
//pub use self::nonlinear_shape_cast_halfspace_support_map::{cast_shapes_nonlinear_halfspace_support_map, cast_shapes_nonlinear_support_map_halfspace};
pub use self::nonlinear_rigid_motion::NonlinearRigidMotion;
pub use self::nonlinear_shape_cast::cast_shapes_nonlinear;
pub use self::nonlinear_shape_cast_support_map_support_map::{
    cast_shapes_nonlinear_support_map_support_map, NonlinearShapeCastMode,
};

#[cfg(feature = "alloc")]
mod nonlinear_shape_cast_composite_shape_shape;
#[cfg(feature = "alloc")]
mod nonlinear_shape_cast_voxels_shape;
//mod cast_shapes_nonlinear_halfspace_support_map;
mod nonlinear_rigid_motion;
mod nonlinear_shape_cast;
mod nonlinear_shape_cast_support_map_support_map;
