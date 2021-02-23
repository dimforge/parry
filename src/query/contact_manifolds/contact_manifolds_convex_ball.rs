use crate::math::{Isometry, Point, Real};
use crate::query::{ContactManifold, TrackedContact};
use crate::shape::{Ball, Shape};
use na::Unit;

/// Computes the contact manifold between a convex shape and a ball, both represented as a `Shape` trait-object.
pub fn contact_manifold_convex_ball_shapes<ManifoldData, ContactData>(
    pos12: &Isometry<Real>,
    shape1: &dyn Shape,
    shape2: &dyn Shape,
    prediction: Real,
    manifold: &mut ContactManifold<ManifoldData, ContactData>,
) where
    ContactData: Default + Copy,
{
    if let Some(ball1) = shape1.as_ball() {
        contact_manifold_convex_ball(&pos12.inverse(), shape2, ball1, prediction, manifold, true);
    } else if let Some(ball2) = shape2.as_ball() {
        contact_manifold_convex_ball(pos12, shape1, ball2, prediction, manifold, false);
    }
}

/// Computes the contact manifold between a convex shape and a ball.
pub fn contact_manifold_convex_ball<'a, ManifoldData, ContactData, S1>(
    pos12: &Isometry<Real>,
    shape1: &'a S1,
    ball2: &'a Ball,
    prediction: Real,
    manifold: &mut ContactManifold<ManifoldData, ContactData>,
    flipped: bool,
) where
    S1: ?Sized + Shape,
    ContactData: Default + Copy,
{
    let local_p2_1 = Point::from(pos12.translation.vector);
    let proj = shape1.project_local_point(&local_p2_1, false);
    let dpos = local_p2_1 - proj.point;

    if let Some((mut local_n1, mut dist)) = Unit::try_new_and_get(dpos, 0.0) {
        if proj.is_inside {
            local_n1 = -local_n1;
            dist = -dist;
        }

        if dist <= ball2.radius + prediction {
            let local_n2 = pos12.inverse_transform_vector(&-*local_n1);
            let local_p2 = (local_n2 * ball2.radius).into();
            let contact_point =
                TrackedContact::flipped(proj.point, local_p2, 0, 0, dist - ball2.radius, flipped);

            if manifold.points.len() != 1 {
                manifold.clear();
                manifold.points.push(contact_point);
            } else {
                // Copy only the geometry so we keep the warmstart impulses.
                manifold.points[0].copy_geometry_from(contact_point);
            }

            if flipped {
                manifold.local_n1 = local_n2;
                manifold.local_n2 = *local_n1;
            } else {
                manifold.local_n1 = *local_n1;
                manifold.local_n2 = local_n2;
            }
        } else {
            manifold.clear();
        }
    }
}
