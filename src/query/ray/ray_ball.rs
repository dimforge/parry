use na::{self, ComplexField};

use crate::math::{Point, Real};
use crate::query::{Ray, RayCast, RayIntersection};
use crate::shape::{Ball, FeatureId};
use num::Zero;

impl RayCast for Ball {
    #[inline]
    fn cast_local_ray(&self, ray: &Ray, max_toi: Real, solid: bool) -> Option<Real> {
        ray_toi_with_ball(&Point::origin(), self.radius, ray, solid)
            .1
            .filter(|toi| *toi <= max_toi)
    }

    #[inline]
    fn cast_local_ray_and_get_normal(
        &self,
        ray: &Ray,
        max_toi: Real,
        solid: bool,
    ) -> Option<RayIntersection> {
        ray_toi_and_normal_with_ball(&Point::origin(), self.radius, ray, solid)
            .1
            .filter(|int| int.toi <= max_toi)
    }
}

/// Computes the time of impact of a ray on a ball.
///
/// The first result element is `true` if the ray started inside of the ball.
#[inline]
pub fn ray_toi_with_ball(
    center: &Point<Real>,
    radius: Real,
    ray: &Ray,
    solid: bool,
) -> (bool, Option<Real>) {
    let dcenter = ray.origin - *center;

    let a = ray.dir.norm_squared();
    let b = dcenter.dot(&ray.dir);
    let c = dcenter.norm_squared() - radius * radius;

    // Special case for when the dir is zero.
    if a.is_zero() {
        if c > 0.0 {
            return (false, None);
        } else {
            return (true, Some(0.0));
        }
    }

    if c > 0.0 && b > 0.0 {
        (false, None)
    } else {
        let delta = b * b - a * c;

        if delta < 0.0 {
            // no solution
            (false, None)
        } else {
            let t = (-b - ComplexField::sqrt(delta)) / a;

            if t <= 0.0 {
                // origin inside of the ball
                if solid {
                    (true, Some(0.0))
                } else {
                    (true, Some((-b + delta.sqrt()) / a))
                }
            } else {
                (false, Some(t))
            }
        }
    }
}

/// Computes the time of impact and contact normal of a ray on a ball.
#[inline]
pub fn ray_toi_and_normal_with_ball(
    center: &Point<Real>,
    radius: Real,
    ray: &Ray,
    solid: bool,
) -> (bool, Option<RayIntersection>) {
    let (inside, inter) = ray_toi_with_ball(&center, radius, ray, solid);

    (
        inside,
        inter.map(|n| {
            let pos = ray.origin + ray.dir * n - center;
            let normal = pos.normalize();

            RayIntersection::new(n, if inside { -normal } else { normal }, FeatureId::Face(0))
        }),
    )
}
