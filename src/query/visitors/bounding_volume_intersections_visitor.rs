use crate::bounding_volume::{SimdAABB, AABB};
use crate::math::SIMD_WIDTH;
use crate::partitioning::{SimdVisitStatus, SimdVisitor};
use simba::simd::SimdBool as _;
use std::marker::PhantomData;

/// Spatial partitioning data structure visitor collecting interferences with a given bounding volume.
pub struct BoundingVolumeIntersectionsVisitor<'a, T: 'a, F> {
    bv: SimdAABB,
    callback: &'a mut F,
    _phantom: PhantomData<T>,
}

impl<'a, T, F> BoundingVolumeIntersectionsVisitor<'a, T, F>
where
    F: FnMut(&T) -> bool,
{
    /// Creates a new `BoundingVolumeIntersectionsVisitor`.
    #[inline]
    pub fn new(bv: &'a AABB, callback: &'a mut F) -> BoundingVolumeIntersectionsVisitor<'a, T, F> {
        BoundingVolumeIntersectionsVisitor {
            bv: SimdAABB::splat(*bv),
            callback,
            _phantom: PhantomData,
        }
    }
}

impl<'a, T, F> SimdVisitor<T, SimdAABB> for BoundingVolumeIntersectionsVisitor<'a, T, F>
where
    F: FnMut(&T) -> bool,
{
    #[inline]
    fn visit(&mut self, bv: &SimdAABB, b: Option<[Option<&T>; SIMD_WIDTH]>) -> SimdVisitStatus {
        let mask = bv.intersects(&self.bv);

        if let Some(data) = b {
            let bitmask = mask.bitmask();

            for ii in 0..SIMD_WIDTH {
                if (bitmask & (1 << ii)) != 0 && data[ii].is_some() {
                    if !(self.callback)(data[ii].unwrap()) {
                        return SimdVisitStatus::ExitEarly;
                    }
                }
            }
        }

        SimdVisitStatus::MaybeContinue(mask)
    }
}
