use crate::bounding_volume::SimdAABB;
use crate::math::{Isometry, Point, Real, SimdBool, SimdReal, Vector, SIMD_WIDTH};
use crate::partitioning::{SimdBestFirstVisitStatus, SimdBestFirstVisitor};
use crate::query::{QueryDispatcher, Ray, SimdRay, TOI};
use crate::shape::{Shape, TypedSimdCompositeShape};
use simba::simd::{SimdBool as _, SimdPartialOrd, SimdValue};

/// Time Of Impact of a composite shape with any other shape, under translational movement.
pub fn time_of_impact_composite_shape_shape<D: ?Sized, G1: ?Sized>(
    dispatcher: &D,
    pos12: &Isometry<Real>,
    vel12: &Vector<Real>,
    g1: &G1,
    g2: &dyn Shape,
    max_toi: Real,
) -> Option<TOI>
where
    D: QueryDispatcher,
    G1: TypedSimdCompositeShape,
{
    let mut visitor =
        TOICompositeShapeShapeBestFirstVisitor::new(dispatcher, pos12, vel12, g1, g2, max_toi);
    g1.typed_qbvh()
        .traverse_best_first(&mut visitor)
        .map(|res| res.1 .1)
}

/// Time Of Impact of any shape with a composite shape, under translational movement.
pub fn time_of_impact_shape_composite_shape<D: ?Sized, G2: ?Sized>(
    dispatcher: &D,
    pos12: &Isometry<Real>,
    vel12: &Vector<Real>,
    g1: &dyn Shape,
    g2: &G2,
    max_toi: Real,
) -> Option<TOI>
where
    D: QueryDispatcher,
    G2: TypedSimdCompositeShape,
{
    time_of_impact_composite_shape_shape(
        dispatcher,
        &pos12.inverse(),
        &-pos12.inverse_transform_vector(&vel12),
        g2,
        g1,
        max_toi,
    )
    .map(|toi| toi.swapped())
}

/// A visitor used to find the time-of-impact between a composite shape and a shape.
pub struct TOICompositeShapeShapeBestFirstVisitor<'a, D: ?Sized, G1: ?Sized + 'a> {
    msum_shift: Vector<SimdReal>,
    msum_margin: Vector<SimdReal>,
    ray: SimdRay,

    dispatcher: &'a D,
    pos12: &'a Isometry<Real>,
    vel12: &'a Vector<Real>,
    g1: &'a G1,
    g2: &'a dyn Shape,
    max_toi: Real,
}

impl<'a, D: ?Sized, G1: ?Sized> TOICompositeShapeShapeBestFirstVisitor<'a, D, G1>
where
    D: QueryDispatcher,
    G1: TypedSimdCompositeShape,
{
    /// Creates a new visitor used to find the time-of-impact between a composite shape and a shape.
    pub fn new(
        dispatcher: &'a D,
        pos12: &'a Isometry<Real>,
        vel12: &'a Vector<Real>,
        g1: &'a G1,
        g2: &'a dyn Shape,
        max_toi: Real,
    ) -> TOICompositeShapeShapeBestFirstVisitor<'a, D, G1> {
        let ls_aabb2 = g2.compute_aabb(pos12);
        let ray = Ray::new(Point::origin(), *vel12);

        TOICompositeShapeShapeBestFirstVisitor {
            dispatcher,
            msum_shift: Vector::splat(-ls_aabb2.center().coords),
            msum_margin: Vector::splat(ls_aabb2.half_extents()),
            ray: SimdRay::splat(ray),
            pos12,
            vel12,
            g1,
            g2,
            max_toi,
        }
    }
}

impl<'a, D: ?Sized, G1: ?Sized> SimdBestFirstVisitor<G1::PartId, SimdAABB>
    for TOICompositeShapeShapeBestFirstVisitor<'a, D, G1>
where
    D: QueryDispatcher,
    G1: TypedSimdCompositeShape,
{
    type Result = (G1::PartId, TOI);

    #[inline]
    fn visit(
        &mut self,
        best: Real,
        bv: &SimdAABB,
        data: Option<[Option<&G1::PartId>; SIMD_WIDTH]>,
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
                    let part_id = *data[ii].unwrap();
                    let mut toi = None;
                    self.g1.map_untyped_part_at(part_id, |part_pos1, g1| {
                        if let Some(part_pos1) = part_pos1 {
                            toi = self
                                .dispatcher
                                .time_of_impact(
                                    &part_pos1.inv_mul(&self.pos12),
                                    &part_pos1.inverse_transform_vector(self.vel12),
                                    g1,
                                    self.g2,
                                    self.max_toi,
                                )
                                .unwrap_or(None)
                                .map(|toi| toi.transform1_by(part_pos1));
                        } else {
                            toi = self
                                .dispatcher
                                .time_of_impact(&self.pos12, self.vel12, g1, self.g2, self.max_toi)
                                .unwrap_or(None);
                        }
                    });

                    if let Some(toi) = toi {
                        results[ii] = Some((part_id, toi));
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
