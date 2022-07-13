use crate::bounding_volume::SimdAABB;
use crate::math::{Isometry, Real, SimdReal, SIMD_WIDTH};
use crate::partitioning::{SimdSimultaneousVisitStatus, SimdSimultaneousVisitor};
use na::SimdValue;
use simba::simd::SimdBool as _;
use std::marker::PhantomData;

/// Spatial partitioning data structure visitor collecting interferences with a given bounding volume.
pub struct BoundingVolumeIntersectionsSimultaneousVisitor<T1, T2, F> {
    pos12: Option<Isometry<SimdReal>>,
    callback: F,
    _phantom: PhantomData<(T1, T2)>,
}

impl<T1, T2, F> BoundingVolumeIntersectionsSimultaneousVisitor<T1, T2, F> {
    /// Creates a new `BoundingVolumeIntersectionsSimultaneousVisitor`.
    #[inline]
    pub fn new(callback: F) -> BoundingVolumeIntersectionsSimultaneousVisitor<T1, T2, F> {
        BoundingVolumeIntersectionsSimultaneousVisitor {
            pos12: None,
            callback,
            _phantom: PhantomData,
        }
    }

    /// Creates a new `BoundingVolumeIntersectionsSimultaneousVisitor`.
    #[inline]
    pub fn with_relative_pos(
        pos12: Isometry<Real>,
        callback: F,
    ) -> BoundingVolumeIntersectionsSimultaneousVisitor<T1, T2, F> {
        BoundingVolumeIntersectionsSimultaneousVisitor {
            pos12: Some(Isometry::splat(pos12)),
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
        let mask = if let Some(pos12) = &self.pos12 {
            let transformed_right_bv = right_bv.transform_by(pos12);
            left_bv.intersects_permutations(&transformed_right_bv)
        } else {
            left_bv.intersects_permutations(right_bv)
        };

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

#[cfg(feature = "parallel")]
impl<T1: Sync, T2: Sync, F> crate::partitioning::ParallelSimdSimultaneousVisitor<T1, T2, SimdAABB>
    for BoundingVolumeIntersectionsSimultaneousVisitor<T1, T2, F>
where
    F: Sync + Fn(&T1, &T2) -> bool,
{
    #[inline]
    fn visit(
        &self,
        left_bv: &SimdAABB,
        left_data: Option<[Option<&T1>; SIMD_WIDTH]>,
        right_bv: &SimdAABB,
        right_data: Option<[Option<&T2>; SIMD_WIDTH]>,
    ) -> SimdSimultaneousVisitStatus {
        let mask = if let Some(pos12) = &self.pos12 {
            let transformed_right_bv = right_bv.transform_by(pos12);
            left_bv.intersects_permutations(&transformed_right_bv)
        } else {
            left_bv.intersects_permutations(right_bv)
        };

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
