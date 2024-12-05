use crate::bounding_volume::{BoundingSphere, SimdAabb};
use crate::math::{Real, SimdBool, SimdReal, SIMD_WIDTH};
use crate::partitioning::{SimdBestFirstVisitStatus, SimdBestFirstVisitor};
use crate::query::{
    self, details::NonlinearShapeCastMode, NonlinearRigidMotion, QueryDispatcher, ShapeCastHit,
};
use crate::shape::{Ball, Shape, TypedSimdCompositeShape};
use simba::simd::SimdValue;

/// Time Of Impact of a composite shape with any other shape, under a rigid motion (translation + rotation).
pub fn cast_shapes_nonlinear_composite_shape_shape<D, G1>(
    dispatcher: &D,
    motion1: &NonlinearRigidMotion,
    g1: &G1,
    motion2: &NonlinearRigidMotion,
    g2: &dyn Shape,
    start_time: Real,
    end_time: Real,
    stop_at_penetration: bool,
) -> Option<ShapeCastHit>
where
    D: ?Sized + QueryDispatcher,
    G1: ?Sized + TypedSimdCompositeShape,
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
pub fn cast_shapes_nonlinear_shape_composite_shape<D, G2>(
    dispatcher: &D,
    motion1: &NonlinearRigidMotion,
    g1: &dyn Shape,
    motion2: &NonlinearRigidMotion,
    g2: &G2,
    start_time: Real,
    end_time: Real,
    stop_at_penetration: bool,
) -> Option<ShapeCastHit>
where
    D: ?Sized + QueryDispatcher,
    G2: ?Sized + TypedSimdCompositeShape,
{
    cast_shapes_nonlinear_composite_shape_shape(
        dispatcher,
        motion2,
        g2,
        motion1,
        g1,
        start_time,
        end_time,
        stop_at_penetration,
    )
    .map(|hit| hit.swapped())
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

impl<'a, D, G1> NonlinearTOICompositeShapeShapeBestFirstVisitor<'a, D, G1>
where
    D: ?Sized + QueryDispatcher,
    G1: ?Sized + TypedSimdCompositeShape,
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

impl<D, G1> SimdBestFirstVisitor<G1::PartId, SimdAabb>
    for NonlinearTOICompositeShapeShapeBestFirstVisitor<'_, D, G1>
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

            if let Some(hit) = query::details::cast_shapes_nonlinear_support_map_support_map(
                self.dispatcher,
                &ball_motion1,
                &ball1,
                &ball1,
                &ball_motion2,
                &ball2,
                &ball2,
                self.start_time,
                self.end_time,
                NonlinearShapeCastMode::StopAtPenetration,
            ) {
                if let Some(data) = data {
                    if hit.time_of_impact < best && data[ii].is_some() {
                        let part_id = *data[ii].unwrap();
                        self.g1.map_untyped_part_at(part_id, |part_pos1, g1, _| {
                            let hit = if let Some(part_pos1) = part_pos1 {
                                self.dispatcher
                                    .cast_shapes_nonlinear(
                                        &self.motion1.prepend(*part_pos1),
                                        g1,
                                        self.motion2,
                                        self.g2,
                                        self.start_time,
                                        self.end_time,
                                        self.stop_at_penetration,
                                    )
                                    .unwrap_or(None)
                                    .map(|hit| hit.transform1_by(part_pos1))
                            } else {
                                self.dispatcher
                                    .cast_shapes_nonlinear(
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

                            // println!("Found time_of_impact: {:?}", time_of_impact);

                            if let Some(hit) = hit {
                                weights[ii] = hit.time_of_impact;
                                mask[ii] = hit.time_of_impact < best;
                                results[ii] = Some((part_id, hit));
                            }
                        });
                    }
                } else {
                    weights[ii] = hit.time_of_impact;
                    mask[ii] = hit.time_of_impact < best;
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
