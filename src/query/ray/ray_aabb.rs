use core::mem;

use crate::bounding_volume::Aabb;
use crate::math::{Real, Vector, VectorExt, DIM};
use crate::query::{Ray, RayCast, RayIntersection};
use crate::shape::FeatureId;
use num::Zero;

impl RayCast for Aabb {
    fn cast_local_ray(&self, ray: &Ray, max_time_of_impact: Real, solid: bool) -> Option<Real> {
        let mut tmin: Real = 0.0;
        let mut tmax: Real = max_time_of_impact;

        for i in 0usize..DIM {
            if ray.dir.vget(i).is_zero() {
                if ray.origin.vget(i) < self.mins.vget(i) || ray.origin.vget(i) > self.maxs.vget(i)
                {
                    return None;
                }
            } else {
                let denom = 1.0 / ray.dir.vget(i);
                let mut inter_with_near_halfspace =
                    (self.mins.vget(i) - ray.origin.vget(i)) * denom;
                let mut inter_with_far_halfspace = (self.maxs.vget(i) - ray.origin.vget(i)) * denom;

                if inter_with_near_halfspace > inter_with_far_halfspace {
                    mem::swap(
                        &mut inter_with_near_halfspace,
                        &mut inter_with_far_halfspace,
                    )
                }

                tmin = tmin.max(inter_with_near_halfspace);
                tmax = tmax.min(inter_with_far_halfspace);

                if tmin > tmax {
                    // This covers the case where tmax is negative because tmin is
                    // initialized at zero.
                    return None;
                }
            }
        }

        if tmin.is_zero() && !solid {
            Some(tmax)
        } else {
            Some(tmin)
        }
    }

    #[inline]
    fn cast_local_ray_and_get_normal(
        &self,
        ray: &Ray,
        max_time_of_impact: Real,
        solid: bool,
    ) -> Option<RayIntersection> {
        ray_aabb(self, ray, max_time_of_impact, solid).map(|(t, n, i)| {
            // NOTE: `i` is the side of the AABB that was hit, as an integer in `[-DIM, DIM]`.
            //       The special value `0` indicates that the ray direction is zero/NaN and its
            //       origin lies inside the AABB; in that case there is no specific face hit.
            let feature = if i == 0 {
                FeatureId::Unknown
            } else if i < 0 {
                FeatureId::Face((-i) as u32 - 1 + 3)
            } else {
                FeatureId::Face(i as u32 - 1)
            };

            RayIntersection::new(t, n, feature)
        })
    }
}

fn ray_aabb(
    aabb: &Aabb,
    ray: &Ray,
    max_time_of_impact: Real,
    solid: bool,
) -> Option<(Real, Vector, isize)> {
    use crate::query::clip;
    clip::clip_aabb_line(aabb, ray.origin, ray.dir).and_then(|(near, far)| {
        if near.0 < 0.0 {
            if solid {
                Some((0.0, Vector::ZERO, far.2))
            } else if far.0 <= max_time_of_impact {
                Some(far)
            } else {
                None
            }
        } else if near.0 <= max_time_of_impact {
            Some(near)
        } else {
            None
        }
    })
}

#[cfg(test)]
mod test {
    use super::*;

    /// Regression test for <https://github.com/dimforge/parry/issues/383>:
    /// casting a ray with a zero direction starting inside an AABB used to panic with
    /// "attempt to subtract with overflow" while computing the hit feature id.
    #[test]
    fn cast_zero_dir_ray_from_inside_aabb_does_not_overflow() {
        let aabb = Aabb::new(-Vector::splat(1.0), Vector::splat(1.0));
        let ray = Ray::new(Vector::ZERO, Vector::ZERO);

        // Both `solid` and non-`solid` variants must not panic.
        let solid = aabb.cast_local_ray_and_get_normal(&ray, 100.0, true);
        let hollow = aabb.cast_local_ray_and_get_normal(&ray, 100.0, false);

        assert!(solid.is_some());
        assert!(hollow.is_some());
        assert_eq!(solid.unwrap().feature, FeatureId::Unknown);
    }
}
