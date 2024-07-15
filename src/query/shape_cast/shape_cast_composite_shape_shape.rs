use crate::bounding_volume::SimdAabb;
use crate::math::{Isometry, Point, Real, SimdBool, SimdReal, Vector, SIMD_WIDTH};
use crate::partitioning::{SimdBestFirstVisitStatus, SimdBestFirstVisitor};
use crate::query::shape_cast::ShapeCastOptions;
use crate::query::{QueryDispatcher, Ray, ShapeCastHit, SimdRay};
use crate::shape::{Shape, TypedSimdCompositeShape};
use simba::simd::{SimdBool as _, SimdPartialOrd, SimdValue};

/// Time Of Impact of a composite shape with any other shape, under translational movement.
pub fn cast_shapes_composite_shape_shape<D, G1>(
    dispatcher: &D,
    pos12: &Isometry<Real>,
    vel12: &Vector<Real>,
    g1: &G1,
    g2: &dyn Shape,
    options: ShapeCastOptions,
) -> Option<ShapeCastHit>
where
    D: ?Sized + QueryDispatcher,
    G1: ?Sized + TypedSimdCompositeShape,
{
    let mut visitor =
        TOICompositeShapeShapeBestFirstVisitor::new(dispatcher, pos12, vel12, g1, g2, options);
    g1.typed_qbvh()
        .traverse_best_first(&mut visitor)
        .map(|res| res.1 .1)
}

/// Time Of Impact of any shape with a composite shape, under translational movement.
pub fn cast_shapes_shape_composite_shape<D, G2>(
    dispatcher: &D,
    pos12: &Isometry<Real>,
    vel12: &Vector<Real>,
    g1: &dyn Shape,
    g2: &G2,
    options: ShapeCastOptions,
) -> Option<ShapeCastHit>
where
    D: ?Sized + QueryDispatcher,
    G2: ?Sized + TypedSimdCompositeShape,
{
    cast_shapes_composite_shape_shape(
        dispatcher,
        &pos12.inverse(),
        &-pos12.inverse_transform_vector(vel12),
        g2,
        g1,
        options,
    )
    .map(|time_of_impact| time_of_impact.swapped())
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
    options: ShapeCastOptions,
}

impl<'a, D, G1> TOICompositeShapeShapeBestFirstVisitor<'a, D, G1>
where
    D: ?Sized + QueryDispatcher,
    G1: ?Sized + TypedSimdCompositeShape,
{
    /// Creates a new visitor used to find the time-of-impact between a composite shape and a shape.
    pub fn new(
        dispatcher: &'a D,
        pos12: &'a Isometry<Real>,
        vel12: &'a Vector<Real>,
        g1: &'a G1,
        g2: &'a dyn Shape,
        options: ShapeCastOptions,
    ) -> TOICompositeShapeShapeBestFirstVisitor<'a, D, G1> {
        let ls_aabb2 = g2.compute_aabb(pos12);
        let ray = Ray::new(Point::origin(), *vel12);

        TOICompositeShapeShapeBestFirstVisitor {
            dispatcher,
            msum_shift: Vector::splat(-ls_aabb2.center().coords),
            msum_margin: Vector::splat(
                ls_aabb2.half_extents() + Vector::repeat(options.target_distance),
            ),
            ray: SimdRay::splat(ray),
            pos12,
            vel12,
            g1,
            g2,
            options,
        }
    }
}

impl<'a, D, G1> SimdBestFirstVisitor<G1::PartId, SimdAabb>
    for TOICompositeShapeShapeBestFirstVisitor<'a, D, G1>
where
    D: ?Sized + QueryDispatcher,
    G1: ?Sized + TypedSimdCompositeShape,
{
    type Result = (G1::PartId, ShapeCastHit);

    #[inline]
    fn visit(
        &mut self,
        best: Real,
        bv: &SimdAabb,
        data: Option<[Option<&G1::PartId>; SIMD_WIDTH]>,
    ) -> SimdBestFirstVisitStatus<Self::Result> {
        // Compute the minkowski sum of the two Aabbs.
        let msum = SimdAabb {
            mins: bv.mins + self.msum_shift + (-self.msum_margin),
            maxs: bv.maxs + self.msum_shift + self.msum_margin,
        };

        // Compute the time of impact.
        let (mask, time_of_impact) =
            msum.cast_local_ray(&self.ray, SimdReal::splat(self.options.max_time_of_impact));

        if let Some(data) = data {
            let better_toi = time_of_impact.simd_lt(SimdReal::splat(best));
            let bitmask = (mask & better_toi).bitmask();
            let mut weights = [0.0; SIMD_WIDTH];
            let mut mask = [false; SIMD_WIDTH];
            let mut results = [None; SIMD_WIDTH];

            for ii in 0..SIMD_WIDTH {
                if (bitmask & (1 << ii)) != 0 && data[ii].is_some() {
                    let part_id = *data[ii].unwrap();
                    let mut hit = None;
                    self.g1.map_untyped_part_at(part_id, |part_pos1, g1, _| {
                        if let Some(part_pos1) = part_pos1 {
                            hit = self
                                .dispatcher
                                .cast_shapes(
                                    &part_pos1.inv_mul(self.pos12),
                                    &part_pos1.inverse_transform_vector(self.vel12),
                                    g1,
                                    self.g2,
                                    self.options,
                                )
                                .unwrap_or(None)
                                .map(|hit| hit.transform1_by(part_pos1));
                        } else {
                            hit = self
                                .dispatcher
                                .cast_shapes(self.pos12, self.vel12, g1, self.g2, self.options)
                                .unwrap_or(None);
                        }
                    });

                    if let Some(hit) = hit {
                        results[ii] = Some((part_id, hit));
                        mask[ii] = hit.time_of_impact < best;
                        weights[ii] = hit.time_of_impact;
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
                weights: time_of_impact,
                mask,
                results: [None; SIMD_WIDTH],
            }
        }
    }
}
