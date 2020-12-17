use crate::math::Real;
use crate::query::Ray;

/// Bounding Volume Tree visitor collecting interferences with a given ray.
pub struct RayInterferencesCollector<'a, T: 'a> {
    /// Ray to be tested.
    pub ray: &'a Ray,
    /// The maximum allowed time of impact.
    pub max_toi: Real,
    /// The data contained by the nodes which bounding volume intersects `self.ray`.
    pub collector: &'a mut Vec<T>,
}

impl<'a, T> RayInterferencesCollector<'a, T> {
    /// Creates a new `RayInterferencesCollector`.
    #[inline]
    pub fn new(
        ray: &'a Ray,
        max_toi: Real,
        buffer: &'a mut Vec<T>,
    ) -> RayInterferencesCollector<'a, T> {
        RayInterferencesCollector {
            ray,
            max_toi,
            collector: buffer,
        }
    }
}

// impl<'a, T, BV> Visitor<T, BV> for RayInterferencesCollector<'a, T>
// where
//     T: Clone,
//     BV: RayCast,
// {
//     #[inline]
//     fn visit(&mut self, bv: &BV, t: Option<&T>) -> VisitStatus {
//         if bv.intersects_local_ray(self.ray, self.max_toi) {
//             if let Some(t) = t {
//                 self.collector.push(t.clone())
//             }
//
//             VisitStatus::Continue
//         } else {
//             VisitStatus::Stop
//         }
//     }
// }
