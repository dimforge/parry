use crate::bounding_volume::{BoundingSphere, SimdAABB};
use crate::math::{Real, SimdBool, SimdReal, SIMD_WIDTH};
use crate::partitioning::{SimdBestFirstVisitStatus, SimdBestFirstVisitor};
use crate::query::{self, details::NonlinearTOIMode, NonlinearRigidMotion, QueryDispatcher, TOI};
use crate::shape::{Ball, Shape, TypedSimdCompositeShape};
use simba::simd::SimdValue;

/// Time Of Impact of a composite shape with any other shape, under a rigid motion (translation + rotation).
pub fn nonlinear_time_of_impact_composite_shape_shape<D: ?Sized, G1: ?Sized>(
    dispatcher: &D,
    motion1: &NonlinearRigidMotion,
    g1: &G1,
    motion2: &NonlinearRigidMotion,
    g2: &dyn Shape,
    start_time: Real,
    end_time: Real,
    stop_at_penetration: bool,
) -> Option<TOI>
where
    D: QueryDispatcher,
    G1: TypedSimdCompositeShape,
{
    let mut visitor = NonlinearTOICompositeShapeShapeBestFirstVisitor::new(
        dispatcher,
        motion1,
        g1,
        motion2,
        g2,
        start_time,
        end_time,
        stop_at_penetration,
    );

    g1.typed_qbvh()
        .traverse_best_first(&mut visitor)
        .map(|res| res.1 .1)
}

/// Time Of Impact of any shape with a composite shape, under a rigid motion (translation + rotation).
pub fn nonlinear_time_of_impact_shape_composite_shape<D: ?Sized, G2: ?Sized>(
    dispatcher: &D,
    motion1: &NonlinearRigidMotion,
    g1: &dyn Shape,
    motion2: &NonlinearRigidMotion,
    g2: &G2,
    start_time: Real,
    end_time: Real,
    stop_at_penetration: bool,
) -> Option<TOI>
where
    D: QueryDispatcher,
    G2: TypedSimdCompositeShape,
{
    nonlinear_time_of_impact_composite_shape_shape(
        dispatcher,
        motion2,
        g2,
        motion1,
        g1,
        start_time,
        end_time,
        stop_at_penetration,
    )
    .map(|toi| toi.swapped())
}

/// A visitor used to determine the non-linear time of impact between a composite shape and another shape.
pub struct NonlinearTOICompositeShapeShapeBestFirstVisitor<'a, D: ?Sized, G1: ?Sized + 'a> {
    sphere2: BoundingSphere,
    start_time: Real,
    end_time: Real,
    stop_at_penetration: bool,

    dispatcher: &'a D,
    motion1: &'a NonlinearRigidMotion,
    motion2: &'a NonlinearRigidMotion,
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
        motion1: &'a NonlinearRigidMotion,
        g1: &'a G1,
        motion2: &'a NonlinearRigidMotion,
        g2: &'a dyn Shape,
        start_time: Real,
        end_time: Real,
        stop_at_penetration: bool,
    ) -> NonlinearTOICompositeShapeShapeBestFirstVisitor<'a, D, G1> {
        NonlinearTOICompositeShapeShapeBestFirstVisitor {
            dispatcher,
            sphere2: g2.compute_local_bounding_sphere(),
            start_time,
            end_time,
            stop_at_penetration,
            motion1,
            motion2,
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
            let ball_motion1 = self.motion1.prepend_translation(center1.coords);
            let ball_motion2 = self.motion2.prepend_translation(self.sphere2.center.coords);

            if let Some(toi) = query::details::nonlinear_time_of_impact_support_map_support_map(
                self.dispatcher,
                &ball_motion1,
                &ball1,
                &ball1,
                &ball_motion2,
                &ball2,
                &ball2,
                self.start_time,
                self.end_time,
                NonlinearTOIMode::StopAtPenetration,
            ) {
                if let Some(data) = data {
                    if toi.toi < best && data[ii].is_some() {
                        let part_id = *data[ii].unwrap();
                        self.g1.map_untyped_part_at(part_id, |part_pos1, g1| {
                            let toi = if let Some(part_pos1) = part_pos1 {
                                self.dispatcher
                                    .nonlinear_time_of_impact(
                                        &self.motion1.prepend(*part_pos1),
                                        g1,
                                        self.motion2,
                                        self.g2,
                                        self.start_time,
                                        self.end_time,
                                        self.stop_at_penetration,
                                    )
                                    .unwrap_or(None)
                                    .map(|toi| toi.transform1_by(part_pos1))
                            } else {
                                self.dispatcher
                                    .nonlinear_time_of_impact(
                                        self.motion1,
                                        g1,
                                        self.motion2,
                                        self.g2,
                                        self.start_time,
                                        self.end_time,
                                        self.stop_at_penetration,
                                    )
                                    .unwrap_or(None)
                            };

                            // println!("Found toi: {:?}", toi);

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
