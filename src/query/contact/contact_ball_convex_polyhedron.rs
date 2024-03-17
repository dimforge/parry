use crate::math::*;
use crate::query::Contact;
use crate::shape::{Ball, Shape};

/// Contact between a ball and a convex polyhedron.
///
/// This function panics if the input shape does not implement
/// both the ConvexPolyhedron and PointQuery traits.
#[inline]
pub fn contact_ball_convex_polyhedron(
    pos12: &Isometry,
    ball1: &Ball,
    shape2: &(impl Shape + ?Sized),
    prediction: Real,
) -> Option<Contact> {
    contact_convex_polyhedron_ball(&pos12.inverse(), shape2, ball1, prediction).map(|c| c.flipped())
}

/// Contact between a convex polyhedron and a ball.
///
/// This function panics if the input shape does not implement
/// both the ConvexPolyhedron and PointQuery traits.
#[inline]
pub fn contact_convex_polyhedron_ball(
    pos12: &Isometry,
    shape1: &(impl Shape + ?Sized),
    ball2: &Ball,
    prediction: Real,
) -> Option<Contact> {
    let center2_1 = Point::from(pos12.translation);
    let (proj, f1) = shape1.project_local_point_and_get_feature(&center2_1);

    let dist;
    let normal1;
    if let Some((dir1, len)) = UnitVector::try_new_and_get(proj.point - center2_1, DEFAULT_EPSILON)
    {
        if proj.is_inside {
            dist = -len - ball2.radius;
            normal1 = dir1;
        } else {
            dist = len - ball2.radius;
            normal1 = -dir1;
        }
    } else {
        dist = -ball2.radius;
        normal1 = shape1
            .feature_normal_at_point(f1, &proj.point)
            .or_else(|| UnitVector::try_new(proj.point.into_vector(), DEFAULT_EPSILON))
            .unwrap_or_else(|| Vector::y_axis());
    }

    if dist <= prediction {
        let normal2 = pos12.inverse_transform_unit_vector(&-normal1);
        let point2 = Point::from(normal2.into_inner() * ball2.radius);
        let point1 = proj.point;
        return Some(Contact::new(point1, point2, normal1, normal2, dist));
    }

    None
}
