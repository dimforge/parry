use crate::math::Real;
use crate::query::{
    DefaultQueryDispatcher, NonlinearRigidMotion, QueryDispatcher, Unsupported, TOI,
};
use crate::shape::Shape;

/// Computes the smallest time of impact of two shapes under arbitrary smooth movement.
pub fn nonlinear_time_of_impact(
    motion1: &NonlinearRigidMotion,
    g1: &dyn Shape,
    motion2: &NonlinearRigidMotion,
    g2: &dyn Shape,
    start_time: Real,
    end_time: Real,
    stop_at_penetration: bool,
) -> Result<Option<TOI>, Unsupported> {
    DefaultQueryDispatcher.nonlinear_time_of_impact(
        motion1,
        g1,
        motion2,
        g2,
        start_time,
        end_time,
        stop_at_penetration,
    )
}
