use crate::bounding_volume::SimdAABB;
use crate::math::{Point, Real};
use crate::partitioning::{SimdVisitStatus, SimdVisitor};
use crate::simd::{SimdReal, SIMD_WIDTH};
use simba::simd::{SimdBool as _, SimdValue};

// FIXME: add a point cost fn.

/// Spatial partitioning structure visitor collecting nodes that may contain a given point.
pub struct PointInterferencesCollector<'a, T: 'a> {
    /// Point to be tested.
    pub point: &'a Point<Real>,
    /// The data contained by the nodes which bounding volume contain `self.point`.
    pub collector: &'a mut Vec<T>,
}

impl<'a, T> PointInterferencesCollector<'a, T> {
    /// Creates a new `PointInterferencesCollector`.
    #[inline]
    pub fn new(
        point: &'a Point<Real>,
        buffer: &'a mut Vec<T>,
    ) -> PointInterferencesCollector<'a, T> {
        PointInterferencesCollector {
            point,
            collector: buffer,
        }
    }
}

impl<'a, T> SimdVisitor<T, SimdAABB> for PointInterferencesCollector<'a, T>
where
    T: Clone,
{
    #[inline]
    fn visit(&mut self, bv: &SimdAABB, b: Option<[Option<&T>; SIMD_WIDTH]>) -> SimdVisitStatus {
        let simd_point: Point<SimdReal> = Point::splat(*self.point);
        let mask = bv.contains_local_point(&simd_point);

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
