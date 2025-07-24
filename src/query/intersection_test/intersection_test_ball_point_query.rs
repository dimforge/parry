use crate::math::{Isometry, Point, Real};
use crate::query::{PointQuery, QueryOptions};
use crate::shape::Ball;

/// Intersection test between a ball and a shape implementing the `PointQuery` trait.
pub fn intersection_test_ball_point_query<P: ?Sized + PointQuery>(
    pos12: &Isometry<Real>,
    ball1: &Ball,
    point_query2: &P,
    options: &dyn QueryOptions,
) -> bool {
    intersection_test_point_query_ball(&pos12.inverse(), point_query2, ball1, options)
}

/// Intersection test between a shape implementing the `PointQuery` trait and a ball.
pub fn intersection_test_point_query_ball<P: ?Sized + PointQuery>(
    pos12: &Isometry<Real>,
    point_query1: &P,
    ball2: &Ball,
    options: &dyn QueryOptions,
) -> bool {
    let local_p2_1 = Point::from(pos12.translation.vector);
    let proj = point_query1.project_local_point(&local_p2_1, true, options);
    proj.is_inside || (local_p2_1 - proj.point).norm_squared() <= ball2.radius * ball2.radius
}
