use crate::bounding_volume::Aabb;
use crate::math::{Point, Real};
use crate::query::{Ray, RayCast, RayIntersection};
use crate::shape::Cuboid;

impl RayCast for Cuboid {
    #[inline]
    fn cast_local_ray(&self, ray: &Ray, max_time_of_impact: Real, solid: bool) -> Option<Real> {
        let dl = Point::from(-self.half_extents);
        let ur = Point::from(self.half_extents);
        Aabb::new(dl, ur).cast_local_ray(ray, max_time_of_impact, solid)
    }

    #[inline]
    fn cast_local_ray_and_get_normal(
        &self,
        ray: &Ray,
        max_time_of_impact: Real,
        solid: bool,
    ) -> Option<RayIntersection> {
        let dl = Point::from(-self.half_extents);
        let ur = Point::from(self.half_extents);
        Aabb::new(dl, ur).cast_local_ray_and_get_normal(ray, max_time_of_impact, solid)
    }
}
