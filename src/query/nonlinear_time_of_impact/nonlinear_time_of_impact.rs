use crate::math::Real;
use crate::motion::RigidMotion;
use crate::query::{DefaultQueryDispatcher, QueryDispatcher, Unsupported, TOI};
use crate::shape::Shape;

/// Computes the smallest time of impact of two shapes under arbitrary smooth movement.
pub fn nonlinear_time_of_impact(
    motion12: &dyn RigidMotion,
    g1: &dyn Shape,
    g2: &dyn Shape,
    max_toi: Real,
    target_distance: Real,
) -> Result<Option<TOI>, Unsupported> {
    DefaultQueryDispatcher.nonlinear_time_of_impact(motion12, g1, g2, max_toi, target_distance)
}
