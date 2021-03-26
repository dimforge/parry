//! Implementation details of the `nonlinear_time_of_impact` function.

pub use self::nonlinear_time_of_impact_composite_shape_shape::{
    nonlinear_time_of_impact_composite_shape_shape, nonlinear_time_of_impact_shape_composite_shape,
    NonlinearTOICompositeShapeShapeBestFirstVisitor,
};
//pub use self::nonlinear_time_of_impact_halfspace_support_map::{nonlinear_time_of_impact_halfspace_support_map, nonlinear_time_of_impact_support_map_halfspace};
pub use self::nonlinear_rigid_motion::NonlinearRigidMotion;
pub use self::nonlinear_time_of_impact::nonlinear_time_of_impact;
pub use self::nonlinear_time_of_impact_support_map_support_map::{
    nonlinear_time_of_impact_support_map_support_map, NonlinearTOIMode,
};

mod nonlinear_time_of_impact_composite_shape_shape;
//mod nonlinear_time_of_impact_halfspace_support_map;
mod nonlinear_rigid_motion;
mod nonlinear_time_of_impact;
mod nonlinear_time_of_impact_support_map_support_map;
