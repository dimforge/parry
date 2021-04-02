use na::Unit;

use crate::math::{Isometry, Point, Real, Vector};
use crate::query::{self, Ray, TOIStatus, TOI};
use crate::shape::Ball;
use num::Zero;

/// Time Of Impact of two balls under translational movement.
#[inline]
pub fn time_of_impact_ball_ball(
    pos12: &Isometry<Real>,
    vel12: &Vector<Real>,
    b1: &Ball,
    b2: &Ball,
    max_toi: Real,
) -> Option<TOI> {
    let rsum = b1.radius + b2.radius;
    let radius = rsum;
    let center = Point::from(-pos12.translation.vector);
    let ray = Ray::new(Point::origin(), *vel12);

    if let (inside, Some(toi)) = query::details::ray_toi_with_ball(&center, radius, &ray, true) {
        if toi > max_toi {
            return None;
        }

        let dpt = ray.point_at(toi) - center;
        let normal1;
        let normal2;
        let witness1;
        let witness2;

        if radius.is_zero() {
            normal1 = Vector::x_axis();
            normal2 = pos12.inverse_transform_unit_vector(&(-Vector::x_axis()));
            witness1 = Point::origin();
            witness2 = Point::origin();
        } else {
            normal1 = Unit::new_unchecked(dpt / radius);
            normal2 = pos12.inverse_transform_unit_vector(&(-normal1));
            witness1 = Point::from(*normal1 * b1.radius);
            witness2 = Point::from(*normal2 * b2.radius);
        }

        let status = if inside && center.coords.norm_squared() < rsum * rsum {
            TOIStatus::Penetrating
        } else {
            TOIStatus::Converged
        };

        Some(TOI {
            toi,
            normal1,
            normal2,
            witness1,
            witness2,
            status,
        })
    } else {
        None
    }
}
