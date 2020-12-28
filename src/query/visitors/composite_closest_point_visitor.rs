use crate::bounding_volume::SimdAABB;
use crate::math::{Point, Real, SimdBool, SimdReal, SIMD_WIDTH};
use crate::partitioning::{SimdBestFirstVisitStatus, SimdBestFirstVisitor};
use crate::query::{PointProjection, PointQuery};
use crate::shape::SimdCompositeShape;
use na;
use simba::simd::{SimdBool as _, SimdPartialOrd, SimdValue};

/// Best-first traversal visitor for computing the point closest to a composite shape.
pub struct CompositeClosestPointVisitor<'a, S: 'a> {
    shape: &'a S,
    point: &'a Point<Real>,
    simd_point: Point<SimdReal>,
    solid: bool,
}

impl<'a, S> CompositeClosestPointVisitor<'a, S> {
    /// Initializes a visitor that allows the computation of the point closest to `point` on `shape`.
    pub fn new(shape: &'a S, point: &'a Point<Real>, solid: bool) -> Self {
        CompositeClosestPointVisitor {
            shape,
            point,
            simd_point: Point::splat(*point),
            solid,
        }
    }
}

impl<'a, S: SimdCompositeShape + PointQuery> SimdBestFirstVisitor<u32, SimdAABB>
    for CompositeClosestPointVisitor<'a, S>
{
    type Result = PointProjection;

    #[inline]
    fn visit(
        &mut self,
        best: Real,
        aabb: &SimdAABB,
        data: Option<[Option<&u32>; SIMD_WIDTH]>,
    ) -> SimdBestFirstVisitStatus<Self::Result> {
        let dist = aabb.distance_to_local_point(&self.simd_point);
        let mask = dist.simd_lt(SimdReal::splat(best));

        if let Some(data) = data {
            let bitmask = mask.bitmask();
            let mut weights = [0.0; SIMD_WIDTH];
            let mut mask = [false; SIMD_WIDTH];
            let mut results = [None; SIMD_WIDTH];

            for ii in 0..SIMD_WIDTH {
                if (bitmask & (1 << ii)) != 0 && data[ii].is_some() {
                    self.shape
                        .map_part_at(*data[ii].unwrap(), &mut |part_pos, obj| {
                            let proj = if let Some(part_pos) = part_pos {
                                obj.project_point(part_pos, self.point, self.solid)
                            } else {
                                obj.project_local_point(self.point, self.solid)
                            };

                            weights[ii] = na::distance(self.point, &proj.point);
                            mask[ii] = true;
                            results[ii] = Some(proj);
                        });
                }
            }

            SimdBestFirstVisitStatus::MaybeContinue {
                weights: SimdReal::from(weights),
                mask: SimdBool::from(mask),
                results,
            }
        } else {
            SimdBestFirstVisitStatus::MaybeContinue {
                weights: dist,
                mask,
                results: [None; SIMD_WIDTH],
            }
        }
    }
}
