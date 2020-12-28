use crate::bounding_volume::{BoundingVolume, SimdAABB};
use crate::math::{Isometry, Point, Real, SimdBool, SimdReal, Vector, SIMD_WIDTH};
use crate::partitioning::{SimdBestFirstVisitStatus, SimdBestFirstVisitor};
use crate::query::{QueryDispatcher, Ray, SimdRay, TOI};
use crate::shape::{Shape, SimdCompositeShape};
use simba::simd::{SimdBool as _, SimdPartialOrd, SimdValue};

/// Time Of Impact of a composite shape with any other shape, under translational movement.
pub fn time_of_impact_composite_shape_shape<D: ?Sized, G1: ?Sized>(
    dispatcher: &D,
    pos12: &Isometry<Real>,
    vel12: &Vector<Real>,
    g1: &G1,
    g2: &dyn Shape,
    max_toi: Real,
    target_distance: Real,
) -> Option<TOI>
where
    D: QueryDispatcher,
    G1: SimdCompositeShape,
{
    let mut visitor = CompositeShapeAgainstAnyTOIVisitor::new(
        dispatcher,
        pos12,
        vel12,
        g1,
        g2,
        max_toi,
        target_distance,
    );
    g1.quadtree()
        .traverse_best_first(&mut visitor)
        .map(|res| res.1)
}

/// Time Of Impact of any shape with a composite shape, under translational movement.
pub fn time_of_impact_shape_composite_shape<D: ?Sized, G2: ?Sized>(
    dispatcher: &D,
    pos12: &Isometry<Real>,
    vel12: &Vector<Real>,
    g1: &dyn Shape,
    g2: &G2,
    max_toi: Real,
    target_distance: Real,
) -> Option<TOI>
where
    D: QueryDispatcher,
    G2: SimdCompositeShape,
{
    time_of_impact_composite_shape_shape(
        dispatcher,
        &pos12.inverse(),
        &-vel12,
        g2,
        g1,
        max_toi,
        target_distance,
    )
    .map(|toi| toi.swapped())
}

struct CompositeShapeAgainstAnyTOIVisitor<'a, D: ?Sized, G1: ?Sized + 'a> {
    msum_shift: Vector<SimdReal>,
    msum_margin: Vector<SimdReal>,
    ray: SimdRay,

    dispatcher: &'a D,
    pos12: &'a Isometry<Real>,
    vel12: &'a Vector<Real>,
    g1: &'a G1,
    g2: &'a dyn Shape,
    max_toi: Real,
    target_distance: Real,
}

impl<'a, D: ?Sized, G1: ?Sized> CompositeShapeAgainstAnyTOIVisitor<'a, D, G1>
where
    D: QueryDispatcher,
    G1: SimdCompositeShape,
{
    pub fn new(
        dispatcher: &'a D,
        pos12: &'a Isometry<Real>,
        vel12: &'a Vector<Real>,
        g1: &'a G1,
        g2: &'a dyn Shape,
        max_toi: Real,
        target_distance: Real,
    ) -> CompositeShapeAgainstAnyTOIVisitor<'a, D, G1> {
        let ls_aabb2 = g2.compute_aabb(pos12).loosened(target_distance);
        let ray = Ray::new(Point::origin(), *vel12);

        CompositeShapeAgainstAnyTOIVisitor {
            dispatcher,
            msum_shift: Vector::splat(-ls_aabb2.center().coords),
            msum_margin: Vector::splat(ls_aabb2.half_extents()),
            ray: SimdRay::splat(ray),
            pos12,
            vel12,
            g1,
            g2,
            max_toi,
            target_distance,
        }
    }
}

impl<'a, D: ?Sized, G1: ?Sized> SimdBestFirstVisitor<u32, SimdAABB>
    for CompositeShapeAgainstAnyTOIVisitor<'a, D, G1>
where
    D: QueryDispatcher,
    G1: SimdCompositeShape,
{
    type Result = TOI;

    #[inline]
    fn visit(
        &mut self,
        best: Real,
        bv: &SimdAABB,
        data: Option<[Option<&u32>; SIMD_WIDTH]>,
    ) -> SimdBestFirstVisitStatus<Self::Result> {
        // Compute the minkowski sum of the two AABBs.
        let msum = SimdAABB {
            mins: bv.mins + self.msum_shift + (-self.msum_margin),
            maxs: bv.maxs + self.msum_shift + self.msum_margin,
        };

        // Compute the TOI.
        let (mask, toi) = msum.cast_local_ray(&self.ray, SimdReal::splat(self.max_toi));

        if let Some(data) = data {
            let better_toi = toi.simd_lt(SimdReal::splat(best));
            let bitmask = (mask & better_toi).bitmask();
            let mut weights = [0.0; SIMD_WIDTH];
            let mut mask = [false; SIMD_WIDTH];
            let mut results = [None; SIMD_WIDTH];

            for ii in 0..SIMD_WIDTH {
                if (bitmask & (1 << ii)) != 0 && data[ii].is_some() {
                    let mut toi = None;
                    self.g1
                        .map_part_at(*data[ii].unwrap(), &mut |part_pos1, g1| {
                            if let Some(part_pos1) = part_pos1 {
                                toi = self
                                    .dispatcher
                                    .time_of_impact(
                                        &part_pos1.inv_mul(&self.pos12),
                                        self.vel12,
                                        g1,
                                        self.g2,
                                        self.max_toi,
                                        self.target_distance,
                                    )
                                    .unwrap_or(None)
                                    .map(|toi| toi.transform1_by(part_pos1));
                            } else {
                                toi = self
                                    .dispatcher
                                    .time_of_impact(
                                        &self.pos12,
                                        self.vel12,
                                        g1,
                                        self.g2,
                                        self.max_toi,
                                        self.target_distance,
                                    )
                                    .unwrap_or(None);
                            }
                        });

                    if let Some(toi) = toi {
                        results[ii] = Some(toi);
                        mask[ii] = toi.toi < best;
                        weights[ii] = toi.toi;
                    }
                }
            }

            SimdBestFirstVisitStatus::MaybeContinue {
                weights: SimdReal::from(weights),
                mask: SimdBool::from(mask),
                results,
            }
        } else {
            SimdBestFirstVisitStatus::MaybeContinue {
                weights: toi,
                mask,
                results: [None; SIMD_WIDTH],
            }
        }
    }
}
