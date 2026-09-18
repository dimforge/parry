use crate::math::Real;
use crate::partitioning::BvhNode;
use crate::query::{
    self, details::NonlinearShapeCastMode, NonlinearRigidMotion, QueryDispatcher, ShapeCastHit,
};
use crate::shape::{Ball, CompositeShapeRef, Shape, SubShapeId, TypedCompositeShape};

impl<S: ?Sized + TypedCompositeShape> CompositeShapeRef<'_, S> {
    /// Performs a non-linear shape-cast between `self` animated subject to the `motion1` and
    /// the `shape2` subject to the `motion2`.
    ///
    /// Returns the index of the sub-shape of `self` that was involved alongside the hit, which is
    /// left as that sub-shape reported it (its `subshape1` is the sub-shape's own when it is a
    /// composite too, and `subshape2` is `shape2`'s).
    pub fn cast_shape_nonlinear<D: ?Sized + QueryDispatcher>(
        &self,
        dispatcher: &D,
        motion1: &NonlinearRigidMotion,
        motion2: &NonlinearRigidMotion,
        shape2: &dyn Shape,
        start_time: Real,
        end_time: Real,
        stop_at_penetration: bool,
    ) -> Option<(SubShapeId, ShapeCastHit)> {
        let sphere2 = shape2.compute_local_bounding_sphere();

        self.0.bvh().find_best(
            end_time,
            |node: &BvhNode, _| {
                let aabb1 = node.aabb();
                let center1 = aabb1.center();
                let radius1 = aabb1.half_extents().length();
                let ball1 = Ball::new(radius1);
                let ball2 = Ball::new(sphere2.radius());
                let ball_motion1 = motion1.prepend_translation(center1);
                let ball_motion2 = motion2.prepend_translation(sphere2.center);

                query::details::cast_shapes_nonlinear_support_map_support_map(
                    dispatcher,
                    &ball_motion1,
                    &ball1,
                    &ball1,
                    &ball_motion2,
                    &ball2,
                    &ball2,
                    start_time,
                    end_time,
                    NonlinearShapeCastMode::StopAtPenetration,
                )
                .map(|hit| hit.time_of_impact)
                .unwrap_or(Real::MAX)
            },
            |part_id, _| {
                self.0
                    .map_untyped_part_at(part_id, |part_pos1, part_shape1, _| {
                        if let Some(part_pos1) = part_pos1 {
                            dispatcher
                                .cast_shapes_nonlinear(
                                    &motion1.prepend(*part_pos1),
                                    part_shape1,
                                    motion2,
                                    shape2,
                                    start_time,
                                    end_time,
                                    stop_at_penetration,
                                )
                                .ok()?
                                .map(|hit| hit.transform1_by(part_pos1))
                        } else {
                            dispatcher
                                .cast_shapes_nonlinear(
                                    motion1,
                                    part_shape1,
                                    motion2,
                                    shape2,
                                    start_time,
                                    end_time,
                                    stop_at_penetration,
                                )
                                .ok()?
                        }
                    })?
            },
        )
    }
}

/// Time Of Impact of a composite shape with any other shape, under a rigid motion (translation + rotation).
pub fn cast_shapes_nonlinear_composite_shape_shape<D, G1>(
    dispatcher: &D,
    motion1: &NonlinearRigidMotion,
    shape1: &G1,
    motion2: &NonlinearRigidMotion,
    shape2: &dyn Shape,
    start_time: Real,
    end_time: Real,
    stop_at_penetration: bool,
) -> Option<ShapeCastHit>
where
    D: ?Sized + QueryDispatcher,
    G1: ?Sized + TypedCompositeShape,
{
    CompositeShapeRef(shape1)
        .cast_shape_nonlinear(
            dispatcher,
            motion1,
            motion2,
            shape2,
            start_time,
            end_time,
            stop_at_penetration,
        )
        .map(|(part_id, mut hit)| {
            hit.subshape1 = part_id;
            hit
        })
}

/// Time Of Impact of any shape with a composite shape, under a rigid motion (translation + rotation).
pub fn cast_shapes_nonlinear_shape_composite_shape<D, G2>(
    dispatcher: &D,
    motion1: &NonlinearRigidMotion,
    shape1: &dyn Shape,
    motion2: &NonlinearRigidMotion,
    shape2: &G2,
    start_time: Real,
    end_time: Real,
    stop_at_penetration: bool,
) -> Option<ShapeCastHit>
where
    D: ?Sized + QueryDispatcher,
    G2: ?Sized + TypedCompositeShape,
{
    cast_shapes_nonlinear_composite_shape_shape(
        dispatcher,
        motion2,
        shape2,
        motion1,
        shape1,
        start_time,
        end_time,
        stop_at_penetration,
    )
    .map(|hit| hit.swapped())
}
