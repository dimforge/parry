use crate::bounding_volume::{SimdAABB, AABB};
use crate::math::SIMD_WIDTH;
use crate::partitioning::{SimdSimultaneousVisitStatus, SimdSimultaneousVisitor};
use simba::simd::SimdBool as _;
use std::marker::PhantomData;

/// Spatial partitioning data structure visitor collecting interferences with a given bounding volume.
pub struct BoundingVolumeIntersectionsSimultaneousVisitor<T1, T2, F> {
    callback: F,
    _phantom: PhantomData<(T1, T2)>,
}

impl<T1, T2, F> BoundingVolumeIntersectionsSimultaneousVisitor<T1, T2, F>
where
    F: FnMut(&T1, &T2) -> bool,
{
    /// Creates a new `BoundingVolumeIntersectionsSimultaneousVisitor`.
    #[inline]
    pub fn new(callback: F) -> BoundingVolumeIntersectionsSimultaneousVisitor<T1, T2, F> {
        BoundingVolumeIntersectionsSimultaneousVisitor {
            callback,
            _phantom: PhantomData,
        }
    }
}

impl<T1, T2, F> SimdSimultaneousVisitor<T1, T2, SimdAABB>
    for BoundingVolumeIntersectionsSimultaneousVisitor<T1, T2, F>
where
    F: FnMut(&T1, &T2) -> bool,
{
    #[inline]
    fn visit(
        &mut self,
        left_bv: &SimdAABB,
        left_data: Option<[Option<&T1>; SIMD_WIDTH]>,
        right_bv: &SimdAABB,
        right_data: Option<[Option<&T2>; SIMD_WIDTH]>,
    ) -> SimdSimultaneousVisitStatus {
        let mask = left_bv.intersects_permutations(right_bv);

        if let (Some(data1), Some(data2)) = (left_data, right_data) {
            for ii in 0..SIMD_WIDTH {
                let bitmask = mask[ii].bitmask();

                for jj in 0..SIMD_WIDTH {
                    if (bitmask & (1 << jj)) != 0 && data1[ii].is_some() && data2[jj].is_some() {
                        if !(self.callback)(data1[ii].unwrap(), data2[jj].unwrap()) {
                            return SimdSimultaneousVisitStatus::ExitEarly;
                        }
                    }
                }
            }
        }

        SimdSimultaneousVisitStatus::MaybeContinue(mask)
    }
}
