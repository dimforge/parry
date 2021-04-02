use crate::math::Real;
use crate::query::{
    DefaultQueryDispatcher, NonlinearRigidMotion, QueryDispatcher, Unsupported, TOI,
};
use crate::shape::Shape;

/// Computes the smallest time of impact of two shapes under translational an rotational movements.
///
/// # Parameters
/// * `motion1` - The motion of the first shape.
/// * `g1` - The first shape involved in the query.
/// * `motion2` - The motion of the second shape.
/// * `g2` - The second shape involved in the query.
/// * `start_time` - The starting time of the interval where the motion takes place.
/// * `end_time` - The end time of the interval where the motion takes place.
/// * `stop_at_penetration` - If the casted shape starts in a penetration state with any
///    collider, two results are possible. If `stop_at_penetration` is `true` then, the
///    result will have a `toi` equal to `start_time`. If `stop_at_penetration` is `false`
///    then the nonlinear shape-casting will see if further motion wrt. the penetration normal
///    would result in tunnelling. If it does not (i.e. we have a separating velocity along
///    that normal) then the nonlinear shape-casting will attempt to find another impact,
///    at a time `> start_time` that could result in tunnelling.
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
