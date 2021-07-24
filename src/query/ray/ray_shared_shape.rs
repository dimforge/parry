use super::{Ray, RayCast, RayIntersection};
use crate::{math::Real, shape::SharedShape};

impl RayCast for SharedShape {
    #[inline]
    fn cast_local_ray(&self, ray: &Ray, max_toi: Real, solid: bool) -> Option<Real> {
        self.0.as_ref().cast_local_ray(ray, max_toi, solid)
    }

    #[inline]
    fn cast_local_ray_and_get_normal(
        &self,
        ray: &Ray,
        max_toi: Real,
        solid: bool,
    ) -> Option<RayIntersection> {
        self.0
            .as_ref()
            .cast_local_ray_and_get_normal(ray, max_toi, solid)
    }
}
