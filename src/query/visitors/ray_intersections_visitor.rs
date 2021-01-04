use crate::bounding_volume::SimdAABB;
use crate::math::{Real, SimdReal, SIMD_WIDTH};
use crate::partitioning::{SimdVisitStatus, SimdVisitor};
use crate::query::{Ray, SimdRay};
use simba::simd::{SimdBool as _, SimdValue};
use std::marker::PhantomData;

/// Bounding Volume Tree visitor collecting intersections with a given ray.
pub struct RayIntersectionsVisitor<'a, T, F> {
    simd_ray: SimdRay,
    max_toi: SimdReal,
    callback: &'a mut F,
    _phantom: PhantomData<T>,
}

impl<'a, T, F> RayIntersectionsVisitor<'a, T, F>
where
    F: FnMut(&T) -> bool,
{
    /// Creates a new `RayIntersectionsVisitor`.
    #[inline]
    pub fn new(ray: &Ray, max_toi: Real, callback: &'a mut F) -> RayIntersectionsVisitor<'a, T, F> {
        RayIntersectionsVisitor {
            simd_ray: SimdRay::splat(*ray),
            max_toi: SimdReal::splat(max_toi),
            callback,
            _phantom: PhantomData,
        }
    }
}

impl<'a, T, F> SimdVisitor<T, SimdAABB> for RayIntersectionsVisitor<'a, T, F>
where
    F: FnMut(&T) -> bool,
{
    #[inline]
    fn visit(&mut self, bv: &SimdAABB, b: Option<[Option<&T>; SIMD_WIDTH]>) -> SimdVisitStatus {
        let mask = bv.cast_local_ray(&self.simd_ray, self.max_toi).0;

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
