use crate::bounding_volume::{BoundingSphere, SimdAABB};
use crate::math::{Real, SimdBool, SimdReal, SIMD_WIDTH};
use crate::motion::{RigidMotion, RigidMotionComposition};
use crate::partitioning::{SimdBestFirstVisitStatus, SimdBestFirstVisitor};
use crate::query::{self, QueryDispatcher, TOI};
use crate::shape::{Ball, Shape, TypedSimdCompositeShape};
use simba::simd::SimdValue;

/// Time Of Impact of a composite shape with any other shape, under a rigid motion (translation + rotation).
pub fn nonlinear_time_of_impact_composite_shape_shape<D: ?Sized, G1: ?Sized>(
    dispatcher: &D,
    motion12: &dyn RigidMotion,
    g1: &G1,
    g2: &dyn Shape,
    max_toi: Real,
    target_distance: Real,
) -> Option<TOI>
where
    D: QueryDispatcher,
    G1: TypedSimdCompositeShape,
{
    let mut visitor = NonlinearTOICompositeShapeShapeBestFirstVisitor::new(
        dispatcher,
        motion12,
        g1,
        g2,
        max_toi,
        target_distance,
    );

    g1.typed_quadtree()
        .traverse_best_first(&mut visitor)
        .map(|res| res.1 .1)
}

/// Time Of Impact of any shape with a composite shape, under a rigid motion (translation + rotation).
pub fn nonlinear_time_of_impact_shape_composite_shape<D: ?Sized, G2: ?Sized>(
    dispatcher: &D,
    motion12: &dyn RigidMotion,
    g1: &dyn Shape,
    g2: &G2,
    max_toi: Real,
    target_distance: Real,
) -> Option<TOI>
where
    D: QueryDispatcher,
    G2: TypedSimdCompositeShape,
{
    nonlinear_time_of_impact_composite_shape_shape(
        dispatcher,
        &motion12.inverse(),
        g2,
        g1,
        max_toi,
        target_distance,
    )
}

/// A visitor used to determine the non-linear time of impact between a composite shape and another shape.
pub struct NonlinearTOICompositeShapeShapeBestFirstVisitor<'a, D: ?Sized, G1: ?Sized + 'a> {
    sphere2: BoundingSphere,
    max_toi: Real,
    target_distance: Real,

    dispatcher: &'a D,
    motion12: &'a dyn RigidMotion,
    g1: &'a G1,
    g2: &'a dyn Shape,
}

impl<'a, D: ?Sized, G1: ?Sized> NonlinearTOICompositeShapeShapeBestFirstVisitor<'a, D, G1>
where
    D: QueryDispatcher,
    G1: TypedSimdCompositeShape,
{
    /// Initializes visitor used to determine the non-linear time of impact between
    /// a composite shape and another shape.
    pub fn new(
        dispatcher: &'a D,
        motion12: &'a dyn RigidMotion,
        g1: &'a G1,
        g2: &'a dyn Shape,
        max_toi: Real,
        target_distance: Real,
    ) -> NonlinearTOICompositeShapeShapeBestFirstVisitor<'a, D, G1> {
        NonlinearTOICompositeShapeShapeBestFirstVisitor {
            dispatcher,
            sphere2: g2.compute_local_aabb().bounding_sphere(),
            max_toi,
            target_distance,
            motion12,
            g1,
            g2,
        }
    }
}

impl<'a, D: ?Sized, G1: ?Sized> SimdBestFirstVisitor<G1::PartId, SimdAABB>
    for NonlinearTOICompositeShapeShapeBestFirstVisitor<'a, D, G1>
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
        let mut weights = [0.0; SIMD_WIDTH];
        let mut mask = [false; SIMD_WIDTH];
        let mut results = [None; SIMD_WIDTH];

        // let centers1: [Point<Real>; SIMD_WIDTH] = bv.center().into();
        let centers1 = bv.center();
        let radius1: [Real; SIMD_WIDTH] = bv.radius().into();

        for ii in 0..SIMD_WIDTH {
            let center1 = centers1.extract(ii);
            let ball1 = Ball::new(radius1[ii]);
            let ball2 = Ball::new(self.sphere2.radius());
            let ball_motion12_part = self
                .motion12
                .prepend_translation(self.sphere2.center().coords);
            let ball_motion12 = ball_motion12_part.append_translation(-center1.coords);

            if let Some(toi) = query::details::nonlinear_time_of_impact_ball_ball(
                &ball_motion12,
                &ball1,
                &ball2,
                self.max_toi,
                self.target_distance,
            ) {
                if let Some(data) = data {
                    if toi.toi < best && data[ii].is_some() {
                        let part_id = *data[ii].unwrap();
                        self.g1.map_untyped_part_at(part_id, |part_pos1, g1| {
                            let toi = if let Some(part_pos1) = part_pos1 {
                                self.dispatcher
                                    .nonlinear_time_of_impact(
                                        &self.motion12.append_transformation(part_pos1.inverse()),
                                        g1,
                                        self.g2,
                                        self.max_toi,
                                        self.target_distance,
                                    )
                                    .unwrap_or(None)
                                    .map(|toi| toi.transform1_by(part_pos1))
                            } else {
                                self.dispatcher
                                    .nonlinear_time_of_impact(
                                        self.motion12,
                                        g1,
                                        self.g2,
                                        self.max_toi,
                                        self.target_distance,
                                    )
                                    .unwrap_or(None)
                            };

                            if let Some(toi) = toi {
                                weights[ii] = toi.toi;
                                mask[ii] = toi.toi < best;
                                results[ii] = Some((part_id, toi));
                            }
                        });
                    }
                } else {
                    weights[ii] = toi.toi;
                    mask[ii] = toi.toi < best;
                }
            }
        }

        SimdBestFirstVisitStatus::MaybeContinue {
            weights: SimdReal::from(weights),
            mask: SimdBool::from(mask),
            results,
        }
    }
}
