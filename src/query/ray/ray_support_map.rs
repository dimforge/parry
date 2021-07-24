use na;

use crate::math::Real;
#[cfg(feature = "dim2")]
use crate::query;
use crate::query::gjk::{self, CSOPoint, VoronoiSimplex};
use crate::query::{Ray, RayCast, RayIntersection};
#[cfg(feature = "dim2")]
use crate::shape::ConvexPolygon;
use crate::shape::{Capsule, FeatureId, Segment, SupportMap};
#[cfg(feature = "dim3")]
use crate::shape::{Cone, ConvexPolyhedron, Cylinder};
use num::Zero;

/// Cast a ray on a shape using the GJK algorithm.
pub fn local_ray_intersection_with_support_map_with_params<G: ?Sized>(
    shape: &G,
    simplex: &mut VoronoiSimplex,
    ray: &Ray,
    max_toi: Real,
    solid: bool,
) -> Option<RayIntersection>
where
    G: SupportMap,
{
    let supp = shape.local_support_point(&-ray.dir);
    simplex.reset(CSOPoint::single_point(supp - ray.origin.coords));

    let inter = gjk::cast_local_ray(shape, simplex, ray, max_toi);

    if !solid {
        inter.and_then(|(toi, normal)| {
            if toi.is_zero() {
                // the ray is inside of the shape.
                let ndir = ray.dir.normalize();
                let supp = shape.local_support_point(&ndir);
                let eps = na::convert::<f64, Real>(0.001f64);
                let shift = (supp - ray.origin).dot(&ndir) + eps;
                let new_ray = Ray::new(ray.origin + ndir * shift, -ray.dir);

                // FIXME: replace by? : simplex.translate_by(&(ray.origin - new_ray.origin));
                simplex.reset(CSOPoint::single_point(supp - new_ray.origin.coords));

                gjk::cast_local_ray(shape, simplex, &new_ray, shift + eps).and_then(
                    |(toi, normal)| {
                        let toi = shift - toi;
                        if toi <= max_toi {
                            Some(RayIntersection::new(toi, normal, FeatureId::Unknown))
                        } else {
                            None
                        }
                    },
                )
            } else {
                Some(RayIntersection::new(toi, normal, FeatureId::Unknown))
            }
        })
    } else {
        inter.map(|(toi, normal)| RayIntersection::new(toi, normal, FeatureId::Unknown))
    }
}

#[cfg(feature = "dim3")]
impl RayCast for Cylinder {
    fn cast_local_ray_and_get_normal(
        &self,
        ray: &Ray,
        max_toi: Real,
        solid: bool,
    ) -> Option<RayIntersection> {
        local_ray_intersection_with_support_map_with_params(
            self,
            &mut VoronoiSimplex::new(),
            &ray,
            max_toi,
            solid,
        )
    }
}

#[cfg(feature = "dim3")]
impl RayCast for Cone {
    fn cast_local_ray_and_get_normal(
        &self,
        ray: &Ray,
        max_toi: Real,
        solid: bool,
    ) -> Option<RayIntersection> {
        local_ray_intersection_with_support_map_with_params(
            self,
            &mut VoronoiSimplex::new(),
            &ray,
            max_toi,
            solid,
        )
    }
}

impl RayCast for Capsule {
    fn cast_local_ray_and_get_normal(
        &self,
        ray: &Ray,
        max_toi: Real,
        solid: bool,
    ) -> Option<RayIntersection> {
        local_ray_intersection_with_support_map_with_params(
            self,
            &mut VoronoiSimplex::new(),
            &ray,
            max_toi,
            solid,
        )
    }
}

#[cfg(feature = "dim3")]
impl RayCast for ConvexPolyhedron {
    fn cast_local_ray_and_get_normal(
        &self,
        ray: &Ray,
        max_toi: Real,
        solid: bool,
    ) -> Option<RayIntersection> {
        local_ray_intersection_with_support_map_with_params(
            self,
            &mut VoronoiSimplex::new(),
            &ray,
            max_toi,
            solid,
        )
    }
}

#[cfg(feature = "dim2")]
impl RayCast for ConvexPolygon {
    fn cast_local_ray_and_get_normal(
        &self,
        ray: &Ray,
        max_toi: Real,
        solid: bool,
    ) -> Option<RayIntersection> {
        local_ray_intersection_with_support_map_with_params(
            self,
            &mut VoronoiSimplex::new(),
            &ray,
            max_toi,
            solid,
        )
    }
}

#[allow(unused_variables)]
impl RayCast for Segment {
    fn cast_local_ray_and_get_normal(
        &self,
        ray: &Ray,
        max_toi: Real,
        solid: bool,
    ) -> Option<RayIntersection> {
        #[cfg(feature = "dim2")]
        {
            use crate::math::Vector;

            let seg_dir = self.scaled_direction();
            let (s, t, parallel) = query::details::closest_points_line_line_parameters_eps(
                &ray.origin,
                &ray.dir,
                &self.a,
                &seg_dir,
                crate::math::DEFAULT_EPSILON,
            );

            if parallel {
                // The lines are parallel, we have to distinguish
                // the case where there is no intersection at all
                // from the case where the line are collinear.
                let dpos = self.a - ray.origin;
                let normal = self.normal().map(|n| *n).unwrap_or_else(Vector::zeros);

                if dpos.dot(&normal).abs() < crate::math::DEFAULT_EPSILON {
                    // The rays and the segment are collinear.
                    let dist1 = dpos.dot(&ray.dir);
                    let dist2 = dist1 + seg_dir.dot(&ray.dir);

                    match (dist1 >= 0.0, dist2 >= 0.0) {
                        (true, true) => {
                            let toi = dist1.min(dist2) / ray.dir.norm_squared();
                            if toi > max_toi {
                                None
                            } else if dist1 <= dist2 {
                                Some(RayIntersection::new(toi, normal, FeatureId::Vertex(0)))
                            } else {
                                Some(RayIntersection::new(
                                    dist2 / ray.dir.norm_squared(),
                                    normal,
                                    FeatureId::Vertex(1),
                                ))
                            }
                        }
                        (true, false) | (false, true) => {
                            // The ray origin lies on the segment.
                            Some(RayIntersection::new(0.0, normal, FeatureId::Face(0)))
                        }
                        (false, false) => {
                            // The segment is behind the ray.
                            None
                        }
                    }
                } else {
                    // The rays never intersect.
                    None
                }
            } else if s >= 0.0 && s <= max_toi && t >= 0.0 && t <= 1.0 {
                let normal = self.normal().map(|n| *n).unwrap_or_else(Vector::zeros);

                if normal.dot(&ray.dir) > 0.0 {
                    Some(RayIntersection::new(s, -normal, FeatureId::Face(1)))
                } else {
                    Some(RayIntersection::new(s, normal, FeatureId::Face(0)))
                }
            } else {
                // The closest points are outside of
                // the ray or segment bounds.
                None
            }
        }
        #[cfg(feature = "dim3")]
        {
            // XXX: implement an analytic solution for 3D too.
            local_ray_intersection_with_support_map_with_params(
                self,
                &mut VoronoiSimplex::new(),
                &ray,
                max_toi,
                solid,
            )
        }
    }
}
