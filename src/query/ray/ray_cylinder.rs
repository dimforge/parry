use crate::math::{ComplexField, Real, Vector};
use crate::query::{Ray, RayCast, RayIntersection};
use crate::shape::{Cylinder, FeatureId};

impl RayCast for Cylinder {
    #[inline]
    fn cast_local_ray_and_get_normal(
        &self,
        ray: &Ray,
        max_time_of_impact: Real,
        solid: bool,
    ) -> Option<RayIntersection> {
        ray_toi_and_normal_with_cylinder(self.half_height, self.radius, ray, solid)
            .filter(|inter| inter.time_of_impact <= max_time_of_impact)
    }
}

/// Computes the intersection of a ray with a cylinder centered at the origin, with its axis
/// along Y.
///
/// The ray is inside the cylinder for the times it is both inside the infinite tube of the
/// given radius and between the two cap planes, so the hit is the later of the two entry
/// times. This is exact, unlike the generic GJK ray cast, which could miss a ray going
/// straight through the middle of a cap when it started far away. Works with an
/// unnormalized ray direction.
pub fn ray_toi_and_normal_with_cylinder(
    half_height: Real,
    radius: Real,
    ray: &Ray,
    solid: bool,
) -> Option<RayIntersection> {
    let (o, d) = (ray.origin, ray.dir);
    if d == Vector::ZERO {
        // Not a ray: it has no direction.
        return None;
    }

    // When the ray is inside the infinite tube: a t² + 2 b t + c <= 0. The quadratic is
    // solved from the point of the ray nearest the axis (`shift` along it), so its terms
    // stay the size of the cylinder: from far away they would nearly cancel.
    let a = d.x * d.x + d.z * d.z;
    let shift = if a == 0.0 {
        0.0
    } else {
        -(o.x * d.x + o.z * d.z) / a
    };
    let near = o + d * shift;
    let b = near.x * d.x + near.z * d.z;
    let c = near.x * near.x + near.z * near.z - radius * radius;
    let (side_in, side_out) = if a == 0.0 {
        // Parallel to the axis: inside the tube for all times, or never.
        if c > 0.0 {
            return None;
        }
        (-Real::MAX, Real::MAX)
    } else {
        let delta = b * b - a * c;
        if delta < 0.0 {
            return None;
        }
        // The numerically stable form of the two roots.
        let root = <Real as ComplexField>::sqrt(delta);
        let q = if b >= 0.0 { -(b + root) } else { root - b };
        if q == 0.0 {
            // b = 0 and delta = 0: the ray only grazes the tube.
            (shift, shift)
        } else {
            let (t0, t1) = (q / a, c / q);
            (shift + t0.min(t1), shift + t0.max(t1))
        }
    };

    // When the ray is between the two cap planes.
    let (cap_in, cap_out) = if d.y == 0.0 {
        if o.y.abs() > half_height {
            return None;
        }
        (-Real::MAX, Real::MAX)
    } else {
        let t0 = (-half_height - o.y) / d.y;
        let t1 = (half_height - o.y) / d.y;
        (t0.min(t1), t0.max(t1))
    };

    let enter = side_in.max(cap_in);
    let exit = side_out.min(cap_out);
    if enter > exit || exit < 0.0 {
        return None;
    }

    // The outward normal where the ray crosses the surface at time `t`, through the cap if
    // the cap plane is the last (entering) or first (leaving) boundary it crosses.
    let normal_at = |t: Real, through_cap: bool, entering: bool| {
        if through_cap {
            let y = if (o.y + d.y * t) >= 0.0 { 1.0 } else { -1.0 };
            Vector::new(0.0, y, 0.0)
        } else {
            // The point on the side, from the nearest point like the roots: `o + d * t` from
            // far away would lose the radius, and could land on the axis itself.
            let p = near + d * (t - shift);
            let across = if entering { -1.0 } else { 1.0 };
            Vector::new(p.x, 0.0, p.z)
                .try_normalize()
                .or_else(|| (Vector::new(d.x, 0.0, d.z) * across).try_normalize())
                .unwrap_or(Vector::X)
        }
    };

    if enter >= 0.0 {
        // The ray starts outside.
        let normal = normal_at(enter, cap_in > side_in, true);
        Some(RayIntersection::new(enter, normal, FeatureId::Unknown))
    } else if solid {
        // The ray starts inside a solid cylinder: as for the other shapes, it hits at once.
        Some(RayIntersection::new(
            0.0,
            normal_at(exit, cap_out < side_out, false),
            FeatureId::Unknown,
        ))
    } else {
        // The ray starts inside a hollow cylinder: it hits the surface on its way out, which
        // faces the ray's origin.
        let normal = -normal_at(exit, cap_out < side_out, false);
        Some(RayIntersection::new(exit, normal, FeatureId::Unknown))
    }
}
