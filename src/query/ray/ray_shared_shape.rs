use crate::{math::Real, shape::SharedShape};

use super::{Ray, RayCast, RayIntersection};


impl RayCast for SharedShape {
  #[inline]
  fn cast_local_ray(&self, ray: &Ray, max_toi: Real, solid: bool) -> Option<Real> {
      if let Some(shape) = self.as_ball() {
        shape.cast_local_ray(ray, max_toi, solid)
      } else if let Some(shape) = self.as_cuboid() {
        shape.cast_local_ray(ray, max_toi, solid)
      } else if let Some(shape) = self.as_halfspace() {
        shape.cast_local_ray(ray, max_toi, solid)
      } else if let Some(shape) = self.as_segment() {
        shape.cast_local_ray(ray, max_toi, solid)
      } else if let Some(shape) = self.as_capsule() {
        shape.cast_local_ray(ray, max_toi, solid)
      } else if let Some(shape) = self.as_triangle() {
        shape.cast_local_ray(ray, max_toi, solid)
      } else if let Some(shape) = self.as_compound() {
        shape.cast_local_ray(ray, max_toi, solid)
      } else if let Some(shape) = self.as_trimesh() {
        shape.cast_local_ray(ray, max_toi, solid)
      } else if let Some(shape) = self.as_polyline() {
        shape.cast_local_ray(ray, max_toi, solid)
      } else if let Some(shape) = self.as_heightfield() {
        shape.cast_local_ray(ray, max_toi, solid)
      } else if let Some(shape) = self.as_round_cuboid() {
        shape.cast_local_ray(ray, max_toi, solid)
      } else if let Some(shape) = self.as_round_triangle() {
        shape.cast_local_ray(ray, max_toi, solid)
      } else {
        #[cfg(feature = "dim3")] {
          if let Some(shape) = self.as_cylinder() {
            shape.cast_local_ray(ray, max_toi, solid)
          } else if let Some(shape) = self.as_cone() {
            shape.cast_local_ray(ray, max_toi, solid)
          } else {
            panic!("Unsupported shape");
          }
        }
        #[cfg(feature = "dim2")] {
          if let Some(shape) = self.as_round_convex_polygon() {
            shape.cast_local_ray(ray, max_toi, solid)
          } else {
            panic!("Unsupported shape");
          }
        }
      }
  }

  #[inline]
  fn cast_local_ray_and_get_normal(
      &self,
      ray: &Ray,
      max_toi: Real,
      solid: bool,
  ) -> Option<RayIntersection> {
    if let Some(shape) = self.as_ball() {
      shape.cast_local_ray_and_get_normal(ray, max_toi, solid)
    } else if let Some(shape) = self.as_cuboid() {
      shape.cast_local_ray_and_get_normal(ray, max_toi, solid)
    } else if let Some(shape) = self.as_halfspace() {
      shape.cast_local_ray_and_get_normal(ray, max_toi, solid)
    } else if let Some(shape) = self.as_segment() {
      shape.cast_local_ray_and_get_normal(ray, max_toi, solid)
    } else if let Some(shape) = self.as_capsule() {
      shape.cast_local_ray_and_get_normal(ray, max_toi, solid)
    } else if let Some(shape) = self.as_triangle() {
      shape.cast_local_ray_and_get_normal(ray, max_toi, solid)
    } else if let Some(shape) = self.as_compound() {
      shape.cast_local_ray_and_get_normal(ray, max_toi, solid)
    } else if let Some(shape) = self.as_trimesh() {
      shape.cast_local_ray_and_get_normal(ray, max_toi, solid)
    } else if let Some(shape) = self.as_polyline() {
      shape.cast_local_ray_and_get_normal(ray, max_toi, solid)
    } else if let Some(shape) = self.as_heightfield() {
      shape.cast_local_ray_and_get_normal(ray, max_toi, solid)
    } else if let Some(shape) = self.as_round_cuboid() {
      shape.cast_local_ray_and_get_normal(ray, max_toi, solid)
    } else if let Some(shape) = self.as_round_triangle() {
      shape.cast_local_ray_and_get_normal(ray, max_toi, solid)
    } else {
      #[cfg(feature = "dim3")] {
        if let Some(shape) = self.as_cylinder() {
          shape.cast_local_ray_and_get_normal(ray, max_toi, solid)
        } else if let Some(shape) = self.as_cone() {
          shape.cast_local_ray_and_get_normal(ray, max_toi, solid)
        } else {
          panic!("Unsupported shape");
        }
      }
      #[cfg(feature = "dim2")] {
        if let Some(shape) = self.as_round_convex_polygon() {
          shape.cast_local_ray_and_get_normal(ray, max_toi, solid)
        } else {
          panic!("Unsupported shape");
        }
      }
    }
  }
}
