use crate::bounding_volume::BoundingSphere;
use crate::math::Real;
use crate::query::{Ray, RayCast, RayIntersection};
use crate::shape::Ball;

impl RayCast for BoundingSphere {
    #[inline]
    fn cast_local_ray(&self, ray: &Ray, max_toi: Real, solid: bool) -> Option<Real> {
        let centered_ray = ray.translate_by(-self.center().coords);
        Ball::new(self.radius()).cast_local_ray(&centered_ray, max_toi, solid)
    }

    #[inline]
    fn cast_local_ray_and_get_normal(
        &self,
        ray: &Ray,
        max_toi: Real,
        solid: bool,
    ) -> Option<RayIntersection> {
        let centered_ray = ray.translate_by(-self.center().coords);
        Ball::new(self.radius()).cast_local_ray_and_get_normal(&centered_ray, max_toi, solid)
    }

    #[inline]
    fn intersects_local_ray(&self, ray: &Ray, max_toi: Real) -> bool {
        let centered_ray = ray.translate_by(-self.center().coords);
        Ball::new(self.radius()).intersects_local_ray(&centered_ray, max_toi)
    }
}
