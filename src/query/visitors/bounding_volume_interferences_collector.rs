use crate::bounding_volume::{SimdAABB, AABB};
use crate::math::SIMD_WIDTH;
use crate::partitioning::{SimdVisitStatus, SimdVisitor};
use simba::simd::SimdBool as _;

/// Spatial partitioning data structure visitor collecting interferences with a given bounding volume.
pub struct BoundingVolumeInterferencesCollector<'a, T: 'a> {
    /// The bounding volume used for interference tests.
    pub bv: SimdAABB,
    /// The data contained by the nodes with bounding volumes intersecting `self.bv`.
    pub collector: &'a mut Vec<T>,
}

impl<'a, T> BoundingVolumeInterferencesCollector<'a, T> {
    /// Creates a new `BoundingVolumeInterferencesCollector`.
    #[inline]
    pub fn new(
        bv: &'a AABB,
        buffer: &'a mut Vec<T>,
    ) -> BoundingVolumeInterferencesCollector<'a, T> {
        BoundingVolumeInterferencesCollector {
            bv: SimdAABB::splat(*bv),
            collector: buffer,
        }
    }
}

impl<'a, T> SimdVisitor<T, SimdAABB> for BoundingVolumeInterferencesCollector<'a, T>
where
    T: Clone,
{
    #[inline]
    fn visit(&mut self, bv: &SimdAABB, b: Option<[Option<&T>; SIMD_WIDTH]>) -> SimdVisitStatus {
        let mask = bv.intersects(&self.bv);

        if let Some(data) = b {
            let bitmask = mask.bitmask();

            for ii in 0..SIMD_WIDTH {
                if (bitmask & (1 << ii)) != 0 && data[ii].is_some() {
                    self.collector.push(data[ii].unwrap().clone())
                }
            }
        }

        SimdVisitStatus::MaybeContinue(mask)
    }
}
