use crate::num::Bounded;
use std::mem;

use na;

use crate::bounding_volume::AABB;
use crate::math::{Point, Real, Vector, DIM};
use crate::query::{Ray, RayCast, RayIntersection};
use crate::shape::{FeatureId, Segment};
use num::Zero;

impl RayCast for AABB {
    fn cast_local_ray(&self, ray: &Ray, max_toi: Real, solid: bool) -> Option<Real> {
        let mut tmin: Real = na::zero::<Real>();
        let mut tmax: Real = max_toi;

        for i in 0usize..DIM {
            if ray.dir[i].is_zero() {
                if ray.origin[i] < self.mins[i] || ray.origin[i] > self.maxs[i] {
                    return None;
                }
            } else {
                let _1: Real = na::one::<Real>();
                let denom = _1 / ray.dir[i];
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

impl AABB {
    /// Computes the parameters of the two intersection points between a line and this AABB.
    ///
    /// The parameters are such that the point are given by `orig + dir * parameter`.
    /// Returns `None` if there is no intersection.
    #[inline]
    pub fn clip_line_parameters(
        &self,
        orig: &Point<Real>,
        dir: &Vector<Real>,
    ) -> Option<(Real, Real)> {
        clip_line(self, orig, dir).map(|clip| ((clip.0).0, (clip.1).0))
    }

    /// Computes the intersection segment between a line and this AABB.
    ///
    /// Returns `None` if there is no intersection.
    #[inline]
    pub fn clip_line(&self, orig: &Point<Real>, dir: &Vector<Real>) -> Option<Segment> {
        clip_line(self, orig, dir)
            .map(|clip| Segment::new(orig + dir * (clip.0).0, orig + dir * (clip.1).0))
    }

    /// Computes the parameters of the two intersection points between a ray and this AABB.
    ///
    /// The parameters are such that the point are given by `ray.orig + ray.dir * parameter`.
    /// Returns `None` if there is no intersection.
    #[inline]
    pub fn clip_ray_parameters(&self, ray: &Ray) -> Option<(Real, Real)> {
        self.clip_line_parameters(&ray.origin, &ray.dir)
            .and_then(|clip| {
                let t0 = clip.0;
                let t1 = clip.1;

                if t1 < na::zero::<Real>() {
                    None
                } else {
                    Some((t0.max(na::zero::<Real>()), t1))
                }
            })
    }

    /// Computes the intersection segment between a ray and this AABB.
    ///
    /// Returns `None` if there is no intersection.
    #[inline]
    pub fn clip_ray(&self, ray: &Ray) -> Option<Segment> {
        self.clip_ray_parameters(ray)
            .map(|clip| Segment::new(ray.point_at(clip.0), ray.point_at(clip.1)))
    }
}

fn clip_line(
    aabb: &AABB,
    origin: &Point<Real>,
    dir: &Vector<Real>,
) -> Option<((Real, Vector<Real>, isize), (Real, Vector<Real>, isize))> {
    // NOTE: we don't start with tmin = 0 so we can return the correct normal
    // when the ray starts exactly on the object contour.

    let mut tmax: Real = Bounded::max_value();
    let mut tmin: Real = -tmax;
    let mut near_side = 0;
    let mut far_side = 0;
    let mut near_diag = false;
    let mut far_diag = false;

    for i in 0usize..DIM {
        if dir[i].is_zero() {
            if origin[i] < aabb.mins[i] || origin[i] > aabb.maxs[i] {
                return None;
            }
        } else {
            let _1: Real = na::one::<Real>();
            let denom = _1 / dir[i];
            let flip_sides;
            let mut inter_with_near_halfspace = (aabb.mins[i] - origin[i]) * denom;
            let mut inter_with_far_halfspace = (aabb.maxs[i] - origin[i]) * denom;

            if inter_with_near_halfspace > inter_with_far_halfspace {
                flip_sides = true;
                mem::swap(
                    &mut inter_with_near_halfspace,
                    &mut inter_with_far_halfspace,
                )
            } else {
                flip_sides = false;
            }

            if inter_with_near_halfspace > tmin {
                tmin = inter_with_near_halfspace;
                near_side = if flip_sides {
                    -(i as isize + 1)
                } else {
                    i as isize + 1
                };
                near_diag = false;
            } else if inter_with_near_halfspace == tmin {
                near_diag = true;
            }

            if inter_with_far_halfspace < tmax {
                tmax = inter_with_far_halfspace;
                far_side = if !flip_sides {
                    -(i as isize + 1)
                } else {
                    i as isize + 1
                };
                far_diag = false;
            } else if inter_with_far_halfspace == tmax {
                far_diag = true;
            }

            if tmax < na::zero::<Real>() || tmin > tmax {
                return None;
            }
        }
    }

    let near = if near_diag {
        (tmin, -dir.normalize(), near_side)
    } else {
        let mut normal = Vector::zeros();

        if near_side < 0 {
            normal[(-near_side - 1) as usize] = na::one::<Real>();
        } else {
            normal[(near_side - 1) as usize] = -na::one::<Real>();
        }

        (tmin, normal, near_side)
    };

    let far = if far_diag {
        (tmax, -dir.normalize(), far_side)
    } else {
        let mut normal = Vector::zeros();

        if far_side < 0 {
            normal[(-far_side - 1) as usize] = -na::one::<Real>();
        } else {
            normal[(far_side - 1) as usize] = na::one::<Real>();
        }

        (tmax, normal, far_side)
    };

    Some((near, far))
}

fn ray_aabb(
    aabb: &AABB,
    ray: &Ray,
    max_toi: Real,
    solid: bool,
) -> Option<(Real, Vector<Real>, isize)> {
    clip_line(aabb, &ray.origin, &ray.dir).and_then(|(near, far)| {
        if near.0 < na::zero::<Real>() {
            if solid {
                Some((na::zero::<Real>(), na::zero(), far.2))
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
