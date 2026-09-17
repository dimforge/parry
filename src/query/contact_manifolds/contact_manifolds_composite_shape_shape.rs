use alloc::{boxed::Box, vec::Vec};

use crate::bounding_volume::BoundingVolume;
use crate::math::{Pose, Real};
use crate::query::contact_manifolds::contact_manifolds_workspace::{
    TypedWorkspaceData, WorkspaceData,
};
use crate::query::contact_manifolds::ContactManifoldsWorkspace;
use crate::query::query_dispatcher::PersistentQueryDispatcher;
use crate::query::ContactManifold;
use crate::shape::{CompositeShape, Shape};
use crate::utils::hashmap::{Entry, HashMap};
use crate::utils::PoseOpt;

#[cfg_attr(feature = "serde-serialize", derive(Serialize, Deserialize))]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Deserialize, rkyv::Serialize)
)]
#[derive(Clone)]
struct SubDetector {
    manifold_id: usize,
    timestamp: bool,
}

#[cfg_attr(feature = "serde-serialize", derive(Serialize, Deserialize))]
#[derive(Clone, Default)]
pub struct CompositeShapeShapeContactManifoldsWorkspace {
    timestamp: bool,
    sub_detectors: HashMap<u32, SubDetector>,
}

impl CompositeShapeShapeContactManifoldsWorkspace {
    pub fn new() -> Self {
        Self::default()
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

/// Leaves reached by the query from which the sub-shape manifolds are computed in parallel
/// (`parallel` feature).
#[cfg(feature = "parallel")]
const PARALLEL_LEAVES: usize = 256;

/// Computes the contact manifolds between a composite shape and an abstract shape.
///
/// The manifolds are in the order of the composite's leaves reached by the query. Under the
/// `parallel` feature, a query reaching many leaves computes their manifolds across threads
/// (same manifolds, same order).
pub fn contact_manifolds_composite_shape_shape<ManifoldData, ContactData>(
    dispatcher: &dyn PersistentQueryDispatcher<ManifoldData, ContactData>,
    pos12: &Pose,
    composite1: &(dyn CompositeShape + Sync),
    shape2: &dyn Shape,
    prediction: Real,
    manifolds: &mut Vec<ContactManifold<ManifoldData, ContactData>>,
    workspace: &mut Option<ContactManifoldsWorkspace>,
    flipped: bool,
) where
    ManifoldData: Default + Clone + Send + Sync,
    ContactData: Default + Copy + Send + Sync,
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
    let deformable = composite1.is_deformable();

    // Traverse bvh1 first.
    let ls_aabb2_1 = shape2.compute_aabb(&pos12).loosened(prediction);
    let mut old_manifolds = core::mem::take(manifolds);

    // The manifold of a leaf: the one kept from the last query, or a fresh one (its sub-shape
    // ids and pose set by the first computation), pushed in leaf order.
    let mut bookkeep = |leaf1: u32, manifolds: &mut Vec<ContactManifold<ManifoldData, ContactData>>| match workspace
        .sub_detectors
        .entry(leaf1)
    {
        Entry::Occupied(entry) => {
            let sub_detector = entry.into_mut();
            let mut manifold = old_manifolds[sub_detector.manifold_id].take();
            sub_detector.manifold_id = manifolds.len();
            sub_detector.timestamp = new_timestamp;
            if deformable {
                manifold.mark_shapes_deformed();
            }
            manifolds.push(manifold);
            false
        }
        Entry::Vacant(entry) => {
            let _ = entry.insert(SubDetector {
                manifold_id: manifolds.len(),
                timestamp: new_timestamp,
            });
            let mut manifold = ContactManifold::new();
            if flipped {
                manifold.subshape1 = 0;
                manifold.subshape2 = leaf1;
            } else {
                manifold.subshape1 = leaf1;
                manifold.subshape2 = 0;
            }
            if deformable {
                manifold.mark_shapes_deformed();
            }
            manifolds.push(manifold);
            true
        }
    };
    // The manifold's sub-shape pose (a fresh manifold) and contacts.
    let compute = |leaf1: u32,
                   fresh: bool,
                   manifold: &mut ContactManifold<ManifoldData, ContactData>| {
        composite1.map_part_at(leaf1, &mut |part_pos1, part_shape1, normal_constraints1| {
            if fresh {
                if flipped {
                    manifold.set_subshape_pos2(part_pos1.copied());
                } else {
                    manifold.set_subshape_pos1(part_pos1.copied());
                }
            }
            if flipped {
                let _ = dispatcher.contact_manifold_convex_convex(
                    &part_pos1.prepend_to(&pos21),
                    shape2,
                    part_shape1,
                    None,
                    normal_constraints1,
                    prediction,
                    manifold,
                );
            } else {
                let _ = dispatcher.contact_manifold_convex_convex(
                    &part_pos1.inv_mul(&pos12),
                    part_shape1,
                    shape2,
                    normal_constraints1,
                    None,
                    prediction,
                    manifold,
                );
            }
        });
    };

    #[cfg(feature = "parallel")]
    {
        let leaves: Vec<u32> = composite1.bvh().intersect_aabb(&ls_aabb2_1).collect();
        if leaves.len() >= PARALLEL_LEAVES {
            use rayon::prelude::*;
            let fresh: Vec<bool> = leaves
                .iter()
                .map(|&leaf1| bookkeep(leaf1, manifolds))
                .collect();
            manifolds
                .par_iter_mut()
                .zip(leaves.par_iter())
                .zip(fresh.par_iter())
                .for_each(|((manifold, &leaf1), &fresh)| compute(leaf1, fresh, manifold));
        } else {
            for leaf1 in leaves {
                let fresh = bookkeep(leaf1, manifolds);
                compute(leaf1, fresh, manifolds.last_mut().unwrap());
            }
        }
    }
    #[cfg(not(feature = "parallel"))]
    for leaf1 in composite1.bvh().intersect_aabb(&ls_aabb2_1) {
        let fresh = bookkeep(leaf1, manifolds);
        compute(leaf1, fresh, manifolds.last_mut().unwrap());
    }

    workspace
        .sub_detectors
        .retain(|_, detector| detector.timestamp == new_timestamp)
}

impl WorkspaceData for CompositeShapeShapeContactManifoldsWorkspace {
    fn as_typed_workspace_data(&self) -> TypedWorkspaceData<'_> {
        TypedWorkspaceData::CompositeShapeShapeContactManifoldsWorkspace(self)
    }

    fn clone_dyn(&self) -> Box<dyn WorkspaceData> {
        Box::new(self.clone())
    }
}
