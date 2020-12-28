use crate::math::Real;
use crate::query::Ray;

/// Bounding Volume Tree visitor collecting intersections with a given ray.
pub struct RayIntersectionsCollector<'a, T: 'a> {
    /// Ray to be tested.
    pub ray: &'a Ray,
    /// The maximum allowed time of impact.
    pub max_toi: Real,
    /// The data contained by the nodes which bounding volume intersects `self.ray`.
    pub collector: &'a mut Vec<T>,
}

impl<'a, T> RayIntersectionsCollector<'a, T> {
    /// Creates a new `RayIntersectionsCollector`.
    #[inline]
    pub fn new(
        ray: &'a Ray,
        max_toi: Real,
        buffer: &'a mut Vec<T>,
    ) -> RayIntersectionsCollector<'a, T> {
        RayIntersectionsCollector {
            ray,
            max_toi,
            collector: buffer,
        }
    }
}

// impl<'a, T, BV> Visitor<T, BV> for RayIntersectionsCollector<'a, T>
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
