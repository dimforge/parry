use crate::math::Real;
use crate::motion::{RigidMotion, RigidMotionComposition};
use crate::query::{DefaultQueryDispatcher, QueryDispatcher, Unsupported, TOI};
use crate::shape::Shape;

/// Computes the smallest time of impact of two shapes under arbitrary smooth movement.
pub fn nonlinear_time_of_impact(
    motion1: &dyn RigidMotion,
    g1: &dyn Shape,
    motion2: &dyn RigidMotion,
    g2: &dyn Shape,
    max_toi: Real,
    target_distance: Real,
) -> Result<Option<TOI>, Unsupported> {
    let motion12 = motion1.inv_mul(motion2);
    DefaultQueryDispatcher.nonlinear_time_of_impact(&motion12, g1, g2, max_toi, target_distance)
}
