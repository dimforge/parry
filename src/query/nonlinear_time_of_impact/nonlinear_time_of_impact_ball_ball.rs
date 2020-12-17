use crate::math::{Isometry, Real};
use crate::motion::RigidMotion;
use crate::query::{self, ClosestPoints, TOI};
use crate::shape::Ball;

/// Non-linear Time Of Impact of two balls under a rigid motion (translation + rotation).
#[inline]
pub fn nonlinear_time_of_impact_ball_ball(
    motion12: &(impl RigidMotion + ?Sized),
    b1: &Ball,
    b2: &Ball,
    max_toi: Real,
    target_distance: Real,
) -> Option<TOI> {
    fn closest_points(
        pos12: &Isometry<Real>,
        g1: &Ball,
        g2: &Ball,
        prediction: Real,
    ) -> ClosestPoints {
        query::details::closest_points_ball_ball(pos12, g1, g2, prediction)
    }

    query::details::nonlinear_time_of_impact_support_map_support_map_with_closest_points_function(
        motion12,
        b1,
        b2,
        max_toi,
        target_distance,
        closest_points,
    )
}
