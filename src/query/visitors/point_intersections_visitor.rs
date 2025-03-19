use crate::bounding_volume::SimdAabb;
use crate::math::{Point, PointT, SimdReal, SIMD_WIDTH};
use crate::partitioning::{SimdVisitStatus, SimdVisitor};
use simba::simd::{SimdBool as _, SimdValue};
use std::marker::PhantomData;

// TODO: add a point cost fn.

/// Spatial partitioning structure visitor collecting nodes that may contain a given point.
pub struct PointIntersectionsVisitor<'a, T, F> {
    simd_point: PointT<SimdReal>,
    /// Callback executed for each leaf which Aabb contains `self.point`.
    callback: &'a mut F,
    _phantom: PhantomData<T>,
}

impl<'a, T, F> PointIntersectionsVisitor<'a, T, F>
where
    F: FnMut(&T) -> bool,
{
    /// Creates a new `PointIntersectionsVisitor`.
    #[inline]
    pub fn new(point: &'a Point, callback: &'a mut F) -> PointIntersectionsVisitor<'a, T, F> {
        PointIntersectionsVisitor {
            simd_point: PointT::splat(*point),
            callback,
            _phantom: PhantomData,
        }
    }
}

impl<T, F> SimdVisitor<T, SimdAabb> for PointIntersectionsVisitor<'_, T, F>
where
    F: FnMut(&T) -> bool,
{
    #[inline]
    fn visit(&mut self, bv: &SimdAabb, b: Option<[Option<&T>; SIMD_WIDTH]>) -> SimdVisitStatus {
        let mask = bv.contains_local_point(&self.simd_point);

        if let Some(data) = b {
            let bitmask = mask.bitmask();

            for (ii, data) in data.iter().enumerate() {
                if (bitmask & (1 << ii)) != 0 {
                    let Some(data) = data else { continue };
                    if !(self.callback)(data) {
                        return SimdVisitStatus::ExitEarly;
                    }
                }
            }
        }

        SimdVisitStatus::MaybeContinue(mask)
    }
}
