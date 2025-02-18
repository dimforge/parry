use crate::math::{Point, SimdReal, Vector, VectorD};
use crate::query::Ray;
use simba::simd::SimdValue;

/// A structure representing 4 rays in an SIMD SoA fashion.
#[derive(Debug, Copy, Clone)]
pub struct SimdRay {
    /// The origin of the rays represented as a single SIMD point.
    pub origin: Point<SimdReal>,
    /// The direction of the rays represented as a single SIMD vector.
    pub dir: VectorD<SimdReal>,
}

impl SimdRay {
    /// Creates a new SIMD ray with all its lanes filled with the same ray.
    pub fn splat(ray: Ray) -> Self {
        Self {
            origin: Point::splat(ray.origin),
            dir: VectorD::splat(ray.dir),
        }
    }
}
