use std::mem;

use na;

use crate::bounding_volume::Aabb;
use crate::math::{Real, Vector, DIM, real};
use crate::query::{Ray, RayCast, RayIntersection};
use crate::shape::FeatureId;
use num::Zero;

impl RayCast for Aabb {
    fn cast_local_ray(&self, ray: &Ray, max_toi: Real, solid: bool) -> Option<Real> {
        let mut tmin: Real = real!(0.0);
        let mut tmax: Real = max_toi;

        for i in 0usize..DIM {
            if ray.dir[i].is_zero() {
                if ray.origin[i] < self.mins[i] || ray.origin[i] > self.maxs[i] {
                    return None;
                }
            } else {
                let denom = real!(1.0) / ray.dir[i];
                let mut inter_with_near_halfspace = (self.mins[i] - ray.origin[i]) * denom;
                let mut inter_with_far_halfspace = (self.maxs[i] - ray.origin[i]) * denom;

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
        max_toi: Real,
        solid: bool,
    ) -> Option<RayIntersection> {
        ray_aabb(self, &ray, max_toi, solid).map(|(t, n, i)| {
            let feature = if i < 0 {
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
    max_toi: Real,
    solid: bool,
) -> Option<(Real, Vector<Real>, isize)> {
    use crate::query::clip;
    clip::clip_aabb_line(aabb, &ray.origin, &ray.dir).and_then(|(near, far)| {
        if near.0 < real!(0.0) {
            if solid {
                Some((real!(0.0), na::zero(), far.2))
            } else if far.0 <= max_toi {
                Some(far)
            } else {
                None
            }
        } else if near.0 <= max_toi {
            Some(near)
        } else {
            None
        }
    })
}
