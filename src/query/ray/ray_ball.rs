use crate::math::{ComplexField, Real, Vector};
use crate::query::{Ray, RayCast, RayIntersection};
use crate::shape::{Ball, FeatureId};
use num::Zero;

impl RayCast for Ball {
    #[inline]
    fn cast_local_ray(&self, ray: &Ray, max_time_of_impact: Real, solid: bool) -> Option<Real> {
        ray_toi_with_ball(Vector::ZERO, self.radius, ray, solid)
            .1
            .filter(|time_of_impact| *time_of_impact <= max_time_of_impact)
    }

    #[inline]
    fn cast_local_ray_and_get_normal(
        &self,
        ray: &Ray,
        max_time_of_impact: Real,
        solid: bool,
    ) -> Option<RayIntersection> {
        ray_toi_and_normal_with_ball(Vector::ZERO, self.radius, ray, solid)
            .1
            .filter(|int| int.time_of_impact <= max_time_of_impact)
    }
}

/// Computes the time of impact of a ray on a ball.
///
/// The first result element is `true` if the ray started inside of the ball.
#[inline]
pub fn ray_toi_with_ball(
    center: Vector,
    radius: Real,
    ray: &Ray,
    solid: bool,
) -> (bool, Option<Real>) {
    let (inside, hit) = ray_toi_and_offset_with_ball(center, radius, ray, solid);
    (inside, hit.map(|(t, _)| t))
}

/// The time of impact, and the hit point relative to the center (formed without cancellation).
///
/// The discriminant is measured at the ray's point of closest approach to the center, as in
/// the capsule and cylinder ray casts: formed from the ray's origin, `|origin - center|² - r²`
/// loses the radius in f32 once the origin is a few thousand radii away.
#[inline]
fn ray_toi_and_offset_with_ball(
    center: Vector,
    radius: Real,
    ray: &Ray,
    solid: bool,
) -> (bool, Option<(Real, Vector)>) {
    let dcenter = ray.origin - center;

    let a = ray.dir.length_squared();
    let b = dcenter.dot(ray.dir);
    let c = dcenter.length_squared() - radius * radius;

    // Special case for when the dir is zero.
    if a.is_zero() {
        if c > 0.0 {
            return (false, None);
        } else {
            return (true, Some((0.0, dcenter)));
        }
    }

    if c > 0.0 && b > 0.0 {
        return (false, None);
    }

    let t_closest = -b / a;
    let closest = dcenter + ray.dir * t_closest;
    let h = radius * radius - closest.length_squared();

    if h < 0.0 && c > 0.0 {
        // no solution
        return (false, None);
    }

    // From inside (c <= 0) the line always crosses the sphere: rounding must not turn it into
    // a miss.
    let half_chord = <Real as ComplexField>::sqrt(h.max(0.0) / a);
    let t = t_closest - half_chord;

    if t <= 0.0 {
        // origin inside of the ball
        if solid {
            (true, Some((0.0, dcenter)))
        } else {
            (
                true,
                Some((t_closest + half_chord, closest + ray.dir * half_chord)),
            )
        }
    } else {
        (false, Some((t, closest - ray.dir * half_chord)))
    }
}

/// Computes the time of impact and contact normal of a ray on a ball.
#[inline]
pub fn ray_toi_and_normal_with_ball(
    center: Vector,
    radius: Real,
    ray: &Ray,
    solid: bool,
) -> (bool, Option<RayIntersection>) {
    let (inside, inter) = ray_toi_and_offset_with_ball(center, radius, ray, solid);

    (
        inside,
        inter.map(|(n, pos)| {
            // `pos` is zero only when the ray starts at the center (solid) or the ball is far
            // below the f32 resolution of the ray: any direction then is as good as another,
            // and facing the ray is the one a hit on a tiny ball would have.
            let normal = pos
                .try_normalize()
                .or_else(|| (-ray.dir).try_normalize())
                .unwrap_or(Vector::ZERO);

            RayIntersection::new(n, if inside { -normal } else { normal }, FeatureId::Face(0))
        }),
    )
}
