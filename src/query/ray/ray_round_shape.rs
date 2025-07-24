use crate::math::Real;
use crate::query::gjk::{GjkOptions, VoronoiSimplex};
use crate::query::{QueryOptions, Ray, RayCast, RayIntersection};
use crate::shape::{RoundShape, SupportMap};

impl<S: SupportMap> RayCast for RoundShape<S> {
    fn cast_local_ray_and_get_normal(
        &self,
        ray: &Ray,
        max_time_of_impact: Real,
        solid: bool,
        options: &dyn QueryOptions,
    ) -> Option<RayIntersection> {
        let options = if let Some(options) = options.as_any().downcast_ref() {
            options
        } else {
            log::warn!("Incorrect option passed to project_local_point: using default options.");
            &GjkOptions::default()
        };
        crate::query::details::local_ray_intersection_with_support_map_with_params(
            self,
            &mut VoronoiSimplex::new(),
            ray,
            max_time_of_impact,
            solid,
            options,
        )
    }
}
