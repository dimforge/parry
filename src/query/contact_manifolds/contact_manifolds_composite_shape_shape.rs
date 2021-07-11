use crate::bounding_volume::BoundingVolume;
use crate::math::{Isometry, Real};
use crate::query::contact_manifolds::contact_manifolds_workspace::{
    TypedWorkspaceData, WorkspaceData,
};
use crate::query::contact_manifolds::ContactManifoldsWorkspace;
use crate::query::query_dispatcher::PersistentQueryDispatcher;
use crate::query::visitors::BoundingVolumeIntersectionsVisitor;
use crate::query::ContactManifold;
use crate::shape::{Shape, SimdCompositeShape};
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
pub struct CompositeShapeShapeContactManifoldsWorkspace {
    timestamp: bool,
    sub_detectors: HashMap<u32, SubDetector>,
}

impl CompositeShapeShapeContactManifoldsWorkspace {
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
            w.0.downcast_ref::<CompositeShapeShapeContactManifoldsWorkspace>()
        })
        .is_some()
    {
        return;
    }

    *workspace = Some(ContactManifoldsWorkspace(Box::new(
        CompositeShapeShapeContactManifoldsWorkspace::new(),
    )));
}

/// Computes the contact manifolds between a composite shape and an abstract shape.
pub fn contact_manifolds_composite_shape_shape<ManifoldData, ContactData>(
    dispatcher: &dyn PersistentQueryDispatcher<ManifoldData, ContactData>,
    pos12: &Isometry<Real>,
    composite1: &dyn SimdCompositeShape,
    shape2: &dyn Shape,
    prediction: Real,
    manifolds: &mut Vec<ContactManifold<ManifoldData, ContactData>>,
    workspace: &mut Option<ContactManifoldsWorkspace>,
    flipped: bool,
) where
    ManifoldData: Default + Clone,
    ContactData: Default + Copy,
{
    ensure_workspace_exists(workspace);
    let workspace: &mut CompositeShapeShapeContactManifoldsWorkspace =
        workspace.as_mut().unwrap().0.downcast_mut().unwrap();
    let new_timestamp = !workspace.timestamp;
    workspace.timestamp = new_timestamp;

    /*
     * Compute interferences.
     */

    let pos12 = *pos12;
    let pos21 = pos12.inverse();

    // Traverse qbvh1 first.
    let ls_aabb2_1 = shape2.compute_aabb(&pos12).loosened(prediction);
    let mut old_manifolds = std::mem::replace(manifolds, Vec::new());

    let mut leaf1_fn = |leaf1: &u32| {
        composite1.map_part_at(*leaf1, &mut |part_pos1, part_shape1| {
            let sub_detector = match workspace.sub_detectors.entry(*leaf1) {
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
                        manifold.subshape1 = 0;
                        manifold.subshape2 = *leaf1;
                        manifold.subshape_pos2 = part_pos1.copied();
                    } else {
                        manifold.subshape1 = *leaf1;
                        manifold.subshape2 = 0;
                        manifold.subshape_pos1 = part_pos1.copied();
                    };

                    manifolds.push(manifold);
                    entry.insert(sub_detector)
                }
            };

            let manifold = &mut manifolds[sub_detector.manifold_id];

            if flipped {
                let _ = dispatcher.contact_manifold_convex_convex(
                    &part_pos1.prepend_to(&pos21),
                    shape2,
                    part_shape1,
                    prediction,
                    manifold,
                );
            } else {
                let _ = dispatcher.contact_manifold_convex_convex(
                    &part_pos1.inv_mul(&pos12),
                    part_shape1,
                    shape2,
                    prediction,
                    manifold,
                );
            }
        });

        true
    };

    let mut visitor1 = BoundingVolumeIntersectionsVisitor::new(&ls_aabb2_1, &mut leaf1_fn);
    composite1.qbvh().traverse_depth_first(&mut visitor1);

    workspace
        .sub_detectors
        .retain(|_, detector| detector.timestamp == new_timestamp)
}

impl WorkspaceData for CompositeShapeShapeContactManifoldsWorkspace {
    fn as_typed_workspace_data(&self) -> TypedWorkspaceData {
        TypedWorkspaceData::CompositeShapeShapeContactManifoldsWorkspace(self)
    }

    fn clone_dyn(&self) -> Box<dyn WorkspaceData> {
        Box::new(self.clone())
    }
}
