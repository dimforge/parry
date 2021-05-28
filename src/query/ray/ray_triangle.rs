use crate::math::Real;
#[cfg(feature = "dim2")]
use crate::math::Vector;
use crate::query::{Ray, RayCast, RayIntersection};
use crate::shape::{FeatureId, Triangle};
#[cfg(feature = "dim3")]
use {crate::math::Point, na::Vector3};

impl RayCast for Triangle {
    #[inline]
    #[cfg(feature = "dim2")]
    fn cast_local_ray_and_get_normal(
        &self,
        ray: &Ray,
        max_toi: Real,
        solid: bool,
    ) -> Option<RayIntersection> {
        let edges = self.edges();

        if solid {
            // Check if ray starts in triangle
            let perp_sign1 = edges[0].scaled_direction().perp(&(ray.origin - edges[0].a)) > 0.0;
            let perp_sign2 = edges[1].scaled_direction().perp(&(ray.origin - edges[1].a)) > 0.0;
            let perp_sign3 = edges[2].scaled_direction().perp(&(ray.origin - edges[2].a)) > 0.0;

            if perp_sign1 == perp_sign2 && perp_sign1 == perp_sign3 {
                return Some(RayIntersection::new(0.0, Vector::y(), FeatureId::Face(0)));
            }
        }

        let mut best = None;
        let mut smallest_toi = Real::MAX;

        for i in 0..3 {
            if let Some(inter) = edges[i].cast_local_ray_and_get_normal(ray, max_toi, solid) {
                if inter.toi < smallest_toi {
                    smallest_toi = inter.toi;
                    best = Some(inter);
                }
            }
        }

        best
    }

    #[inline]
    #[cfg(feature = "dim3")]
    fn cast_local_ray_and_get_normal(
        &self,
        ray: &Ray,
        max_toi: Real,
        _: bool,
    ) -> Option<RayIntersection> {
        let inter = local_ray_intersection_with_triangle(&self.a, &self.b, &self.c, &ray)?.0;

        if inter.toi <= max_toi {
            Some(inter)
        } else {
            None
        }
    }
}

/// Computes the intersection between a triangle and a ray.
///
/// If an intersection is found, the time of impact, the normal and the barycentric coordinates of
/// the intersection point are returned.
#[cfg(feature = "dim3")]
pub fn local_ray_intersection_with_triangle(
    a: &Point<Real>,
    b: &Point<Real>,
    c: &Point<Real>,
    ray: &Ray,
) -> Option<(RayIntersection, Vector3<Real>)> {
    let ab = *b - *a;
    let ac = *c - *a;

    // normal
    let n = ab.cross(&ac);
    let d = n.dot(&ray.dir);

    // the normal and the ray direction are parallel
    if d == 0.0 {
        return None;
    }

    let ap = ray.origin - *a;
    let t = ap.dot(&n);

    // the ray does not intersect the halfspace defined by the triangle
    if (t < 0.0 && d < 0.0) || (t > 0.0 && d > 0.0) {
        return None;
    }

    let fid = if d < 0.0 { 0 } else { 1 };

    let d = d.abs();

    //
    // intersection: compute barycentric coordinates
    //
    let e = -ray.dir.cross(&ap);

    let mut v;
    let mut w;
    let toi;
    let normal;

    if t < 0.0 {
        v = -ac.dot(&e);

        if v < 0.0 || v > d {
            return None;
        }

        w = ab.dot(&e);

        if w < 0.0 || v + w > d {
            return None;
        }

        let invd = 1.0 / d;
        toi = -t * invd;
        normal = -n.normalize();
        v = v * invd;
        w = w * invd;
    } else {
        v = ac.dot(&e);

        if v < 0.0 || v > d {
            return None;
        }

        w = -ab.dot(&e);

        if w < 0.0 || v + w > d {
            return None;
        }

        let invd = 1.0 / d;
        toi = t * invd;
        normal = n.normalize();
        v = v * invd;
        w = w * invd;
    }

    Some((
        RayIntersection::new(toi, normal, FeatureId::Face(fid)),
        Vector3::new(-v - w + 1.0, v, w),
    ))
}
