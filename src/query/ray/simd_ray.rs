use crate::math::{Point, SimdReal, Vector};
use crate::query::Ray;
use simba::simd::SimdValue;

#[derive(Debug, Copy, Clone)]
pub struct SimdRay {
    pub origin: Point<SimdReal>,
    pub dir: Vector<SimdReal>,
}

impl SimdRay {
    pub fn splat(ray: Ray) -> Self {
        Self {
            origin: Point::splat(ray.origin),
            dir: Vector::splat(ray.dir),
        }
    }
}
