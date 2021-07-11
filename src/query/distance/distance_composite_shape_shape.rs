use crate::bounding_volume::SimdAABB;
use crate::math::{Isometry, Real, SimdBool, SimdReal, Vector, SIMD_WIDTH};
use crate::partitioning::{SimdBestFirstVisitStatus, SimdBestFirstVisitor};
use crate::query::QueryDispatcher;
use crate::shape::{Shape, TypedSimdCompositeShape};
use crate::utils::IsometryOpt;
use simba::simd::{SimdBool as _, SimdPartialOrd, SimdValue};

/// Smallest distance between a composite shape and any other shape.
pub fn distance_composite_shape_shape<D: ?Sized, G1: ?Sized>(
    dispatcher: &D,
    pos12: &Isometry<Real>,
    g1: &G1,
    g2: &dyn Shape,
) -> Real
where
    D: QueryDispatcher,
    G1: TypedSimdCompositeShape,
{
    let mut visitor = CompositeShapeAgainstAnyDistanceVisitor::new(dispatcher, pos12, g1, g2);
    g1.typed_qbvh()
        .traverse_best_first(&mut visitor)
        .expect("The composite shape must not be empty.")
        .1
         .1
}

/// Smallest distance between a shape and a composite shape.
pub fn distance_shape_composite_shape<D: ?Sized, G2: ?Sized>(
    dispatcher: &D,
    pos12: &Isometry<Real>,
    g1: &dyn Shape,
    g2: &G2,
) -> Real
where
    D: QueryDispatcher,
    G2: TypedSimdCompositeShape,
{
    distance_composite_shape_shape(dispatcher, &pos12.inverse(), g2, g1)
}

/// A visitor for computing the distance between a composite shape and a shape.
pub struct CompositeShapeAgainstAnyDistanceVisitor<'a, D: ?Sized, G1: ?Sized + 'a> {
    msum_shift: Vector<SimdReal>,
    msum_margin: Vector<SimdReal>,

    dispatcher: &'a D,
    pos12: &'a Isometry<Real>,
    g1: &'a G1,
    g2: &'a dyn Shape,
}

impl<'a, D: ?Sized, G1: ?Sized + 'a> CompositeShapeAgainstAnyDistanceVisitor<'a, D, G1> {
    /// Initialize a visitor for computing the distance between a composite shape and a shape.
    pub fn new(
        dispatcher: &'a D,
        pos12: &'a Isometry<Real>,
        g1: &'a G1,
        g2: &'a dyn Shape,
    ) -> Self {
        let ls_aabb2 = g2.compute_aabb(pos12);

        Self {
            dispatcher,
            msum_shift: Vector::splat(-ls_aabb2.center().coords),
            msum_margin: Vector::splat(ls_aabb2.half_extents()),
            pos12,
            g1,
            g2,
        }
    }
}

impl<'a, D: ?Sized, G1: ?Sized> SimdBestFirstVisitor<G1::PartId, SimdAABB>
    for CompositeShapeAgainstAnyDistanceVisitor<'a, D, G1>
where
    D: QueryDispatcher,
    G1: TypedSimdCompositeShape,
{
    type Result = (G1::PartId, Real);

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
        let dist = msum.distance_to_origin();
        let mask = dist.simd_lt(SimdReal::splat(best));

        if let Some(data) = data {
            let bitmask = mask.bitmask();
            let mut weights = [0.0; SIMD_WIDTH];
            let mut mask = [false; SIMD_WIDTH];
            let mut results = [None; SIMD_WIDTH];

            for ii in 0..SIMD_WIDTH {
                if (bitmask & (1 << ii)) != 0 && data[ii].is_some() {
                    let part_id = *data[ii].unwrap();
                    let mut dist = Ok(0.0);
                    self.g1.map_untyped_part_at(part_id, |part_pos1, g1| {
                        dist =
                            self.dispatcher
                                .distance(&part_pos1.inv_mul(self.pos12), g1, self.g2);
                    });

                    match dist {
                        Ok(dist) => {
                            if dist == 0.0 {
                                return SimdBestFirstVisitStatus::ExitEarly(Some((part_id, 0.0)));
                            } else {
                                weights[ii] = dist;
                                mask[ii] = dist < best;
                                results[ii] = Some((part_id, dist));
                            }
                        }
                        Err(_) => {}
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
                weights: dist,
                mask,
                results: [None; SIMD_WIDTH],
            }
        }
    }
}
