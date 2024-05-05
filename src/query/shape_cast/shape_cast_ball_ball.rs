use na::Unit;

use crate::math::{Isometry, Point, Real, Vector};
use crate::query::details::ShapeCastOptions;
use crate::query::{self, Ray, ShapeCastHit, ShapeCastStatus};
use crate::shape::Ball;
use num::Zero;

/// Time Of Impact of two balls under translational movement.
#[inline]
pub fn cast_shapes_ball_ball(
    pos12: &Isometry<Real>,
    vel12: &Vector<Real>,
    b1: &Ball,
    b2: &Ball,
    options: ShapeCastOptions,
) -> Option<ShapeCastHit> {
    let rsum = b1.radius + b2.radius + options.target_distance;
    let radius = rsum;
    let center = Point::from(-pos12.translation.vector);
    let ray = Ray::new(Point::origin(), *vel12);

    if let (inside, Some(time_of_impact)) =
        query::details::ray_toi_with_ball(&center, radius, &ray, true)
    {
        if time_of_impact > options.max_time_of_impact {
            return None;
        }

        let dpt = ray.point_at(time_of_impact) - center;
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
            ShapeCastStatus::PenetratingOrWithinTargetDist
        } else {
            ShapeCastStatus::Converged
        };

        Some(ShapeCastHit {
            time_of_impact,
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
