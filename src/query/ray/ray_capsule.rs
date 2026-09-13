use crate::math::{ComplexField, Real, Vector};
use crate::query::{Ray, RayCast, RayIntersection};
use crate::shape::{Capsule, FeatureId, Segment};

impl RayCast for Capsule {
    #[inline]
    fn cast_local_ray(&self, ray: &Ray, max_time_of_impact: Real, solid: bool) -> Option<Real> {
        ray_toi_with_capsule(&self.segment, self.radius, ray, solid)
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
        ray_toi_and_normal_with_capsule(&self.segment, self.radius, ray, solid)
            .filter(|inter| inter.time_of_impact <= max_time_of_impact)
    }
}

/// Computes the time of impact of a ray on a capsule.
/// Returns true if the ray started inside the capsule and the time of impact.
///
/// Adapted from Inigo Quilez (<https://iquilezles.org/articles/intersectors/>).
/// Adapted to unnormalized ray direction.
/// Made robust to degenerate cases and ray origin inside the capsule.
/// Switched to projecting onto the plane with cross products
/// because they introduce much less error than a difference of dot products.
fn ray_toi_with_capsule(
    segment: &Segment,
    radius: Real,
    ray: &Ray,
    solid: bool,
) -> (bool, Option<Real>) {
    let ab = segment.b - segment.a;
    let ao = ray.origin - segment.a;

    let ab_ab = ab.length_squared();
    let dir_dir = ray.dir.length_squared();
    let ab_dir = ab.dot(ray.dir);
    let ab_ao = ab.dot(ao);
    let radius_squared = radius * radius;

    // do a circle intersection on the plane perpendicular to the capsule's axis.
    // all these variables are scaled by ab
    let dir_on_plane = cross(ray.dir, ab);
    let origin_on_plane = cross(ao, ab);
    let ray_step = dir_on_plane.length_squared();
    let b = dir_on_plane.dot(origin_on_plane);
    let separation = origin_on_plane.length_squared() - radius_squared * ab_ab;
    let h = b * b - ray_step * separation;

    let inside = separation <= 0.0
        && (0.0 < ab_ao || ao.length_squared() <= radius_squared)
        && (ab_ao < ab_ab || (ray.origin - segment.b).length_squared() <= radius_squared);

    if inside && solid {
        return (true, Some(0.0));
    }

    if h >= 0.0 {
        let check_sphere_a = if ray_step == 0.0 {
            // the ray is parallel to the capsule,
            // so it can only hit one of the caps
            (ab_dir > 0.0) ^ inside
        } else {
            // cylinder part
            // when outside, take the first intersection, when inside, the second
            let radical = <Real as ComplexField>::sqrt(h);
            let t = (-b + if inside { radical } else { -radical }) / ray_step;
            let y = ab_ao + t * ab_dir;
            if 0.0 < y && y < ab_ab && t >= 0.0 {
                return (inside, Some(t));
            }
            y <= 0.0
        };

        // caps
        let oc = if check_sphere_a {
            ao
        } else {
            ray.origin - segment.b
        };
        let b = ray.dir.dot(oc);
        let c = oc.length_squared() - radius_squared;
        let h = b * b - c * dir_dir;
        if h >= 0.0 {
            let radical = <Real as ComplexField>::sqrt(h);
            let t = -b + if inside { radical } else { -radical };
            if t >= 0.0 && dir_dir != 0.0 {
                return (inside, Some(t / dir_dir));
            }
        }
    }
    (inside, None)
}

#[cfg(feature = "dim3")]
#[inline]
fn cross(v: Vector, segment: Vector) -> Vector {
    v.cross(segment)
}

/// Returns a vector with zero y, which is complete nonsense
/// but makes the 2D case work with the same code as the 3D case.
#[cfg(feature = "dim2")]
#[inline]
fn cross(v: Vector, segment: Vector) -> Vector {
    Vector::new(v.x * segment.y - v.y * segment.x, 0.0)
}

/// Computes the time of impact and contact normal of a ray on a capsule.
fn ray_toi_and_normal_with_capsule(
    segment: &Segment,
    radius: Real,
    ray: &Ray,
    solid: bool,
) -> Option<RayIntersection> {
    let (inside, inter) = ray_toi_with_capsule(segment, radius, ray, solid);

    inter.map(|t| {
        let normal = if solid && inside {
            Vector::ZERO
        } else {
            let p = ray.origin + ray.dir * t;
            let a_to_p = p - segment.a;
            let seg = segment.b - segment.a;
            let seg_squared = seg.length_squared();

            // the projection of the point onto the capsule's axis times the segment's length
            let proj_times_seg = a_to_p.dot(seg);

            let n = if proj_times_seg <= 0.0 {
                (a_to_p).normalize()
            } else if proj_times_seg >= seg_squared {
                (p - segment.b).normalize()
            } else {
                (a_to_p - (proj_times_seg / seg_squared) * seg).normalize()
            };
            if inside {
                -n
            } else {
                n
            }
        };
        RayIntersection::new(t, normal, FeatureId::Face(0))
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::Vector;
    use crate::query::point::point_query::PointQuery;
    use oorandom::Rand32;

    #[test]
    fn exact_cases() {
        let c = Capsule::new(v2(0.0, 0.5), v2(0.0, 1.5), 0.5);
        // Hit straight down the axis on the top cap, unnormalized direction.
        expect_hit(&c, v2(0.0, 5.0), v2(0.0, -0.2), true, 15.0, v2(0.0, 1.0));
        // Oblique hit on the cylinder.
        expect_hit(&c, v2(5.0, 1.0), v2(-0.3, 0.0), true, 15.0, v2(1.0, 0.0));
        // Tangential hit at the tip of the top cap.
        expect_hit(&c, v2(5.0, 2.0), v2(-1.0, 0.0), true, 5.0, v2(0.0, 1.0));
        // Hit on the bottom cap from below, parallel to the axis but offset
        // from it.
        expect_hit(
            &c,
            v2(0.1, -4.0),
            v2(0.0, 0.2),
            true,
            20.0505,
            v2(0.2, -0.9798),
        );
        // Lateral miss.
        assert!(c
            .cast_local_ray(&Ray::new(v2(10.0, 5.0), v2(0.0, 0.1)), 50.0, true)
            .is_none());

        // Outside-origin misses where the ray's line crosses the tube (or a cap
        // sphere) only BEHIND the origin: the entry root is negative and must
        // not be reported. The fuzz can't catch these (it only aims rays at
        // interior points, so its entries are always positive).
        // Inside the infinite tube past the b-cap, receding. The tube entry is
        // behind (t = -0.9) with a phantom axis coordinate inside the band.
        assert!(c
            .cast_local_ray(&Ray::new(v2(0.4, 1.85), v2(1.0, 1.0)), 50.0, true)
            .is_none());
        // On-axis past the b-cap, parallel, receding (phantom t = -2.1).
        assert!(c
            .cast_local_ray(&Ray::new(v2(0.0, 2.1), v2(0.0, 1.0)), 50.0, true)
            .is_none());
        // Inside the tube past the b-cap, receding at an angle. Phantom axis
        // coordinate beyond the slab, cap entry behind (t = -0.87).
        assert!(c
            .cast_local_ray(&Ray::new(v2(0.4, 1.85), v2(1.0, 0.2)), 50.0, true)
            .is_none());
        // Outside everything, receding. The line crosses the tube behind the
        // origin, phantom axis coordinate in the band (t = -2.5).
        assert!(c
            .cast_local_ray(&Ray::new(v2(2.0, 1.0), v2(1.0, 0.0)), 50.0, true)
            .is_none());
        // Same, with the phantom axis coordinate exactly on the a-end boundary.
        assert!(c
            .cast_local_ray(&Ray::new(v2(2.0, 3.0), v2(1.0, 1.0)), 50.0, true)
            .is_none());
        // Inside, solid: contact at the origin, zero normal.
        expect_hit(&c, v2(0.1, 1.0), v2(0.0, 1.0), true, 0.0, v2(0.0, 0.0));
        // Inside, hollow: the exit, inward normal.
        expect_hit(&c, v2(0.0, 1.0), v2(0.0, 1.0), false, 1.0, v2(0.0, -1.0));
        // Same, toward the a-cap (the other parallel routing branch).
        expect_hit(&c, v2(0.0, 1.0), v2(0.0, -1.0), false, 1.0, v2(0.0, 1.0));
        // Inside, hollow: exit through the cylinder's side (the band exit
        // root, not a cap).
        expect_hit(&c, v2(0.0, 1.0), v2(1.0, 0.0), false, 0.5, v2(-1.0, 0.0));
        // Same, oblique.
        expect_hit(&c, v2(0.1, 0.8), v2(1.0, 0.5), false, 0.4, v2(-1.0, 0.0));
        // Inside the b-cap sphere past the slab: the b-sphere exit (t = 0.6)
        // is an intermediate crossing and must be skipped in favor of the
        // last one (the a-sphere exit at t = 1.6).
        expect_hit(&c, v2(0.0, 1.6), v2(0.0, -1.0), false, 1.6, v2(0.0, 1.0));
        // Degenerate zero-length ray, inside / outside.
        expect_hit(&c, v2(0.1, 1.0), v2(0.0, 0.0), true, 0.0, v2(0.0, 0.0));
        assert!(c
            .cast_local_ray(&Ray::new(v2(0.1, 3.0), v2(0.0, 0.0)), 50.0, true)
            .is_none());
        // Degenerate capsule (a == b): behaves as a ball of radius 1 at (0, 1).
        let ball = Capsule::new(v2(0.0, 1.0), v2(0.0, 1.0), 1.0);
        expect_hit(&ball, v2(0.0, 5.0), v2(0.0, -1.0), true, 3.0, v2(0.0, 1.0));
        expect_hit(&ball, v2(0.5, 1.0), v2(1.0, 0.0), true, 0.0, v2(0.0, 0.0));
        // max_toi filtering (the top-cap hit above is at t = 15).
        assert!(c
            .cast_local_ray(&Ray::new(v2(0.0, 5.0), v2(0.0, -0.2)), 14.9, true)
            .is_none());
        assert!(c
            .cast_local_ray(&Ray::new(v2(0.0, 5.0), v2(0.0, -0.2)), 15.1, true)
            .is_some());
        assert!(c
            .cast_local_ray_and_get_normal(&Ray::new(v2(0.0, 5.0), v2(0.0, -0.2)), 14.9, true)
            .is_none());
    }

    fn v2(x: Real, y: Real) -> Vector {
        Vector::new(
            x,
            y,
            #[cfg(feature = "dim3")]
            0.0,
        )
    }

    #[track_caller]
    fn expect_hit(c: &Capsule, o: Vector, d: Vector, solid: bool, et: Real, en: Vector) {
        let i = c
            .cast_local_ray_and_get_normal(&Ray::new(o, d), 50.0, solid)
            .unwrap_or_else(|| panic!("expected hit (o={:?}, d={:?})", o, d));
        assert!(
            (i.time_of_impact - et).abs() < 1e-4,
            "t: got {}, want {}",
            i.time_of_impact,
            et
        );
        assert!(
            (i.normal - en).length() < 1e-3,
            "n: got {:?}, want {:?}",
            i.normal,
            en
        );
    }

    #[test]
    fn fuzz_capsule_ray_casts() {
        let epsilon = 0.002;
        let mut rng = Rand32::new(42);

        for _ in 0..100_000 {
            let (a, b) = (rnd_vec(&mut rng, 10.0), rnd_vec(&mut rng, 10.0));
            let r = 0.5 + 5.0 * rnd(&mut rng);
            let capsule = Capsule::new(a, b, r);

            // a random point inside the capsule
            let inside = {
                let mut w = rnd_vec(&mut rng, r);
                while w.length_squared() >= r * r {
                    w = rnd_vec(&mut rng, r);
                }
                a + (b - a) * rnd(&mut rng) + w
            };

            // cast random ray toward the inside point
            let far_enough = r + a.distance(b);
            let mut offset = Vector::ZERO;
            while offset.length_squared() < far_enough * far_enough {
                offset = rnd_vec(&mut rng, far_enough * 2.0);
            }

            let o = (a + b) * 0.5 + offset;
            let d = (inside - o) * (0.1 + 0.9 * rnd(&mut rng));
            let i = capsule
                .cast_local_ray_and_get_normal(&Ray::new(o, d), 1000.0, true)
                .expect("a ray aimed at an interior point must hit");

            let hit = o + d * i.time_of_impact;
            assert!(
                capsule.contains_local_point(hit - i.normal * epsilon),
                "nudging inward along the normal should go inside the capsule"
            );
            assert!(
                !capsule.contains_local_point(hit + i.normal * epsilon),
                "nudging outward along the normal should go outside the capsule"
            );

            // A ray from the interior point to the far point (hollow) must
            // exit the capsule; its normal points inward.
            let i_in = capsule
                .cast_local_ray_and_get_normal(&Ray::new(inside, o - inside), 1000.0, false)
                .expect("a ray from inside toward the outside must exit");
            let hit_in = inside + (o - inside) * i_in.time_of_impact;
            assert!(
                capsule.contains_local_point(hit_in + i_in.normal * epsilon),
                "nudging along the inward normal should stay inside"
            );
            assert!(
                !capsule.contains_local_point(hit_in - i_in.normal * epsilon),
                "nudging against the inward normal should go outside"
            );

            assert!(
                capsule
                    .cast_local_ray(&Ray::new(o, -d), 1000.0, true)
                    .is_none(),
                "a retreating ray must miss"
            );

            #[cfg(feature = "dim2")]
            let tangent = Vector::new(-i.normal.y, i.normal.x);
            #[cfg(feature = "dim3")]
            let tangent = {
                let mut tangent = Vector::ZERO;
                while tangent.length_squared() < 1e-8 {
                    tangent = rnd_vec(&mut rng, 1.0).cross(i.normal);
                }
                tangent
            };
            let origin = hit + i.normal * (epsilon + rnd(&mut rng)) - rnd(&mut rng) * tangent;
            assert!(
                capsule
                    .cast_local_ray(&Ray::new(origin, tangent), 1000.0, true)
                    .is_none(),
                "tangent outside the capsule should miss"
            );
        }
    }

    fn rnd(rng: &mut Rand32) -> Real {
        #[cfg(feature = "f32")]
        {
            rng.rand_float()
        }
        #[cfg(feature = "f64")]
        {
            rng.rand_float() as Real
        }
    }

    fn rnd_vec(rng: &mut Rand32, scale: Real) -> Vector {
        let mut component = || (rnd(rng) - 0.5) * 2.0 * scale;
        #[cfg(feature = "dim2")]
        {
            Vector::new(component(), component())
        }
        #[cfg(feature = "dim3")]
        {
            Vector::new(component(), component(), component())
        }
    }
}
