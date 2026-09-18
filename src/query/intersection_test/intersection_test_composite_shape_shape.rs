use crate::bounding_volume::BoundingVolume;
use crate::math::Pose;
use crate::partitioning::BvhNode;
use crate::query::{QueryDispatcher, ShapeIntersection};
use crate::shape::{CompositeShapeRef, Shape, SubShapeId, TypedCompositeShape};
use crate::utils::PoseOpt;

impl<S: ?Sized + TypedCompositeShape> CompositeShapeRef<'_, S> {
    /// Tests whether the given other `shape`, positioned at `pose12` relative to `self`,
    /// intersects `self`.
    ///
    /// Returns `None` when they do not intersect. Otherwise returns the index of a sub-shape of
    /// `self` that `shape` intersects alongside that sub-shape's own result (its `subshape1` is
    /// the sub-shape's own when it is a composite too, and `subshape2` is `shape`'s).
    pub fn intersects_shape<D: ?Sized + QueryDispatcher>(
        &self,
        dispatcher: &D,
        pose12: &Pose,
        shape: &dyn Shape,
    ) -> Option<(SubShapeId, ShapeIntersection)> {
        let ls_aabb2 = shape.compute_aabb(pose12);
        self.0
            .bvh()
            .leaves(|node: &BvhNode| node.aabb().intersects(&ls_aabb2))
            .find_map(|leaf_id| {
                self.0
                    .map_untyped_part_at(leaf_id, |part_pose1, sub1, _| {
                        dispatcher
                            .intersection_test(&part_pose1.inv_mul(pose12), sub1, shape)
                            .ok()
                            .filter(|result| result.intersecting)
                            .map(|result| (leaf_id, result))
                    })
                    .flatten()
            })
    }
}

/// Intersection test between a composite shape (`Mesh`, `Compound`) and any other shape.
pub fn intersection_test_composite_shape_shape<D, G1>(
    dispatcher: &D,
    pos12: &Pose,
    g1: &G1,
    g2: &dyn Shape,
) -> ShapeIntersection
where
    D: ?Sized + QueryDispatcher,
    G1: ?Sized + TypedCompositeShape,
{
    // `subshape2` is left as the dispatch set it: `g2` may be a composite too, and only it
    // knows which of its parts answered.
    match CompositeShapeRef(g1).intersects_shape(dispatcher, pos12, g2) {
        Some((part_id, result)) => {
            ShapeIntersection::new(true).with_subshapes(part_id, result.subshape2)
        }
        None => ShapeIntersection::new(false),
    }
}

/// Proximity between a shape and a composite (`Mesh`, `Compound`) shape.
pub fn intersection_test_shape_composite_shape<D, G2>(
    dispatcher: &D,
    pos12: &Pose,
    g1: &dyn Shape,
    g2: &G2,
) -> ShapeIntersection
where
    D: ?Sized + QueryDispatcher,
    G2: ?Sized + TypedCompositeShape,
{
    intersection_test_composite_shape_shape(dispatcher, &pos12.inverse(), g2, g1).swapped()
}
