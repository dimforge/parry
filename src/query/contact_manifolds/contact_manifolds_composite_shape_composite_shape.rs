use crate::bounding_volume::BoundingVolume;
use crate::math::{Isometry, Real};
use crate::query::contact_manifolds::contact_manifolds_workspace::{
    TypedWorkspaceData, WorkspaceData,
};
use crate::query::contact_manifolds::ContactManifoldsWorkspace;
use crate::query::query_dispatcher::PersistentQueryDispatcher;
use crate::query::visitors::BoundingVolumeIntersectionsVisitor;
use crate::query::ContactManifold;
use crate::shape::SimdCompositeShape;
use crate::utils::hashmap::{Entry, HashMap};
use crate::utils::IsometryOpt;

#[cfg_attr(feature = "serde-serialize", derive(Serialize, Deserialize))]
#[derive(Clone)]
struct SubDetector {
    manifold_id: usize,
    timestamp: bool,
}

#[cfg_attr(feature = "serde-serialize", derive(Serialize, Deserialize))]
#[derive(Clone)]
pub struct CompositeShapeCompositeShapeContactManifoldsWorkspace {
    timestamp: bool,
    sub_detectors: HashMap<(u32, u32), SubDetector>,
}

impl CompositeShapeCompositeShapeContactManifoldsWorkspace {
    pub fn new() -> Self {
        Self {
            timestamp: false,
            sub_detectors: HashMap::default(),
        }
    }
}

fn ensure_workspace_exists(workspace: &mut Option<ContactManifoldsWorkspace>) {
    if workspace
        .as_ref()
        .and_then(|w| {
            w.0.downcast_ref::<CompositeShapeCompositeShapeContactManifoldsWorkspace>()
        })
        .is_some()
    {
        return;
    }

    *workspace = Some(ContactManifoldsWorkspace(Box::new(
        CompositeShapeCompositeShapeContactManifoldsWorkspace::new(),
    )));
}

/// Computes the contact manifolds between two composite shapes.
pub fn contact_manifolds_composite_shape_composite_shape<'a, ManifoldData, ContactData>(
    dispatcher: &dyn PersistentQueryDispatcher<ManifoldData, ContactData>,
    pos12: &Isometry<Real>,
    mut composite1: &'a dyn SimdCompositeShape,
    mut composite2: &'a dyn SimdCompositeShape,
    prediction: Real,
    manifolds: &mut Vec<ContactManifold<ManifoldData, ContactData>>,
    workspace: &mut Option<ContactManifoldsWorkspace>,
) where
    ManifoldData: Default + Clone,
    ContactData: Default + Copy,
{
    ensure_workspace_exists(workspace);
    let workspace: &mut CompositeShapeCompositeShapeContactManifoldsWorkspace =
        workspace.as_mut().unwrap().0.downcast_mut().unwrap();
    let new_timestamp = !workspace.timestamp;
    workspace.timestamp = new_timestamp;

    /*
     * Compute interferences.
     */

    let mut qbvh1 = composite1.qbvh();
    let mut qbvh2 = composite2.qbvh();

    let mut pos12 = *pos12;
    let mut pos21 = pos12.inverse();
    let mut stack2 = Vec::new();

    let mut ls_aabb1 = qbvh1.root_aabb();
    let mut ls_aabb2 = qbvh2.root_aabb();
    let flipped = ls_aabb1.half_extents().norm_squared() < ls_aabb2.half_extents().norm_squared();

    if flipped {
        std::mem::swap(&mut composite1, &mut composite2);
        std::mem::swap(&mut qbvh1, &mut qbvh2);
        std::mem::swap(&mut pos12, &mut pos21);
        std::mem::swap(&mut ls_aabb1, &mut ls_aabb2);
    }

    // Traverse qbvh1 first.
    let ls_aabb2_1 = ls_aabb2.transform_by(&pos12).loosened(prediction);
    let mut old_manifolds = std::mem::replace(manifolds, Vec::new());

    let mut leaf_fn1 = |leaf1: &u32| {
        composite1.map_part_at(*leaf1, &mut |part_pos1, part_shape1| {
            let pos211 = part_pos1.prepend_to(&pos21); // == pos21 * part_pos1
            let ls_part_aabb1_2 = part_shape1.compute_aabb(&pos211).loosened(prediction);
            let mut leaf_fn2 = |leaf2: &u32| {
                composite2.map_part_at(*leaf2, &mut |part_pos2, part_shape2| {
                    let pos2211 = part_pos2.inv_mul(&pos211);
                    let entry_key = if flipped {
                        (*leaf2, *leaf1)
                    } else {
                        (*leaf1, *leaf2)
                    };

                    let sub_detector = match workspace.sub_detectors.entry(entry_key) {
                        Entry::Occupied(entry) => {
                            let sub_detector = entry.into_mut();
                            let manifold = old_manifolds[sub_detector.manifold_id].take();
                            sub_detector.manifold_id = manifolds.len();
                            sub_detector.timestamp = new_timestamp;
                            manifolds.push(manifold);
                            sub_detector
                        }
                        Entry::Vacant(entry) => {
                            let sub_detector = SubDetector {
                                manifold_id: manifolds.len(),
                                timestamp: new_timestamp,
                            };

                            let mut manifold = ContactManifold::new();

                            if flipped {
                                manifold.subshape1 = *leaf2;
                                manifold.subshape2 = *leaf1;
                                manifold.subshape_pos1 = part_pos2.copied();
                                manifold.subshape_pos2 = part_pos1.copied();
                            } else {
                                manifold.subshape1 = *leaf1;
                                manifold.subshape2 = *leaf2;
                                manifold.subshape_pos1 = part_pos1.copied();
                                manifold.subshape_pos2 = part_pos2.copied();
                            };

                            manifolds.push(manifold);
                            entry.insert(sub_detector)
                        }
                    };

                    let manifold = &mut manifolds[sub_detector.manifold_id];

                    if flipped {
                        let _ = dispatcher.contact_manifold_convex_convex(
                            &pos2211,
                            part_shape2,
                            part_shape1,
                            prediction,
                            manifold,
                        );
                    } else {
                        let _ = dispatcher.contact_manifold_convex_convex(
                            &pos2211.inverse(),
                            part_shape1,
                            part_shape2,
                            prediction,
                            manifold,
                        );
                    }
                });

                true
            };

            let mut visitor2 =
                BoundingVolumeIntersectionsVisitor::new(&ls_part_aabb1_2, &mut leaf_fn2);

            qbvh2.traverse_depth_first_with_stack(&mut visitor2, &mut stack2);
        });

        true
    };

    let mut visitor1 = BoundingVolumeIntersectionsVisitor::new(&ls_aabb2_1, &mut leaf_fn1);
    qbvh1.traverse_depth_first(&mut visitor1);

    workspace
        .sub_detectors
        .retain(|_, detector| detector.timestamp == new_timestamp)
}

impl WorkspaceData for CompositeShapeCompositeShapeContactManifoldsWorkspace {
    fn as_typed_workspace_data(&self) -> TypedWorkspaceData {
        TypedWorkspaceData::CompositeShapeCompositeShapeContactManifoldsWorkspace(self)
    }

    fn clone_dyn(&self) -> Box<dyn WorkspaceData> {
        Box::new(self.clone())
    }
}
