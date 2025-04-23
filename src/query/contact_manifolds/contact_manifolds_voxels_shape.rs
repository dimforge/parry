use crate::bounding_volume::Aabb;
use crate::math::{Isometry, Point, Real, Translation, Vector, DIM};
use crate::query::{
    ContactManifold, ContactManifoldsWorkspace, PersistentQueryDispatcher, PointQuery,
    TypedWorkspaceData, WorkspaceData,
};
use crate::shape::{
    AxisMask, Cuboid, RoundShape, Shape, SupportMap, VoxelPrimitiveGeometry, VoxelType, Voxels,
};
use crate::utils::hashmap::{Entry, HashMap};
use crate::utils::IsometryOpt;
use alloc::{boxed::Box, vec::Vec};

#[cfg_attr(feature = "serde-serialize", derive(Serialize, Deserialize))]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Deserialize, rkyv::Serialize),
    archive(check_bytes)
)]
#[derive(Clone)]
struct SubDetector {
    manifold_id: usize,
    selected_contacts: u32,
    timestamp: bool,
}

// NOTE: this is using a similar kind of cache as compound shape and height-field.
//       It is different from the trimesh cash though. Which one is better?
/// A workspace for collision-detection against voxels shape.
#[cfg_attr(feature = "serde-serialize", derive(Serialize, Deserialize))]
#[derive(Clone, Default)]
pub struct VoxelsShapeContactManifoldsWorkspace {
    timestamp: bool,
    sub_detectors: HashMap<[u32; 2], SubDetector>,
}

impl VoxelsShapeContactManifoldsWorkspace {
    /// A new empty workspace for collision-detection against voxels shape.
    pub fn new() -> Self {
        Self::default()
    }
}

impl WorkspaceData for VoxelsShapeContactManifoldsWorkspace {
    fn as_typed_workspace_data(&self) -> TypedWorkspaceData {
        TypedWorkspaceData::VoxelsShapeContactManifoldsWorkspace(self)
    }

    fn clone_dyn(&self) -> Box<dyn WorkspaceData> {
        Box::new(self.clone())
    }
}

fn ensure_workspace_exists(workspace: &mut Option<ContactManifoldsWorkspace>) {
    if workspace
        .as_ref()
        .and_then(|w| w.0.downcast_ref::<VoxelsShapeContactManifoldsWorkspace>())
        .is_some()
    {
        return;
    }

    *workspace = Some(ContactManifoldsWorkspace(Box::new(
        VoxelsShapeContactManifoldsWorkspace::new(),
    )));
}

/// Computes the contact manifold between a convex shape and a ball, both represented as a `Shape` trait-object.
pub fn contact_manifolds_voxels_shape_shapes<ManifoldData, ContactData>(
    dispatcher: &dyn PersistentQueryDispatcher<ManifoldData, ContactData>,
    pos12: &Isometry<Real>,
    shape1: &dyn Shape,
    shape2: &dyn Shape,
    prediction: Real,
    manifolds: &mut Vec<ContactManifold<ManifoldData, ContactData>>,
    workspace: &mut Option<ContactManifoldsWorkspace>,
) where
    ManifoldData: Default + Clone,
    ContactData: Default + Copy,
{
    if let Some(voxels1) = shape1.as_voxels() {
        contact_manifolds_voxels_shape(
            dispatcher, pos12, voxels1, shape2, prediction, manifolds, workspace, false,
        );
    } else if let Some(voxels2) = shape2.as_voxels() {
        contact_manifolds_voxels_shape(
            dispatcher,
            &pos12.inverse(),
            voxels2,
            shape1,
            prediction,
            manifolds,
            workspace,
            true,
        );
    }
}

/// Computes the contact manifold between a convex shape and a ball.
pub fn contact_manifolds_voxels_shape<ManifoldData, ContactData>(
    dispatcher: &dyn PersistentQueryDispatcher<ManifoldData, ContactData>,
    pos12: &Isometry<Real>,
    voxels1: &Voxels,
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
    let workspace: &mut VoxelsShapeContactManifoldsWorkspace =
        workspace.as_mut().unwrap().0.downcast_mut().unwrap();
    let new_timestamp = !workspace.timestamp;

    // TODO: avoid reallocating the new `manifolds` vec at each step.
    let mut old_manifolds = core::mem::take(manifolds);

    workspace.timestamp = new_timestamp;

    let radius1 = voxels1.voxel_size() / 2.0;

    let aabb1 = voxels1.local_aabb();
    let aabb2_1 = shape2.compute_aabb(pos12);
    let domain2_1 = Aabb {
        mins: aabb2_1.mins - radius1 * 10.0,
        maxs: aabb2_1.maxs + radius1 * 10.0,
    };

    if let Some(intersection_aabb1) = aabb1.intersection(&aabb2_1) {
        for (vid, center1, data1) in voxels1.voxels_intersecting_local_aabb(&intersection_aabb1) {
            let voxel1 = data1.voxel_type();

            // TODO: would be nice to have a strategy to handle interior voxels for depenetration.
            if voxel1 == VoxelType::Empty || voxel1 == VoxelType::Interior {
                continue;
            }

            let key_voxel = voxels1.voxel_key_at(vid);
            let mut key_low = key_voxel;
            let mut key_high = key_low;
            let mins = voxels1.domain_mins;
            let maxs = voxels1.domain_maxs - Vector::repeat(1);
            let mask1 = data1.free_faces();

            let adjust_canon = |axis: AxisMask, i: usize, key: &mut Point<i32>, val: i32| {
                if !mask1.contains(axis) {
                    key[i] = val;
                }
            };

            adjust_canon(AxisMask::X_POS, 0, &mut key_high, maxs[0]);
            adjust_canon(AxisMask::X_NEG, 0, &mut key_low, mins[0]);
            adjust_canon(AxisMask::Y_POS, 1, &mut key_high, maxs[1]);
            adjust_canon(AxisMask::Y_NEG, 1, &mut key_low, mins[1]);

            #[cfg(feature = "dim3")]
            {
                adjust_canon(AxisMask::Z_POS, 2, &mut key_high, maxs[2]);
                adjust_canon(AxisMask::Z_NEG, 2, &mut key_low, mins[2]);
            }

            let workspace_key = [
                voxels1.linear_index(key_low),
                voxels1.linear_index(key_high),
            ];

            // TODO: could we refactor the workspace system between Voxels, HeightField, and CompoundShape?
            //       (and maybe TriMesh too but it’s using a different approach).
            let (sub_detector, manifold_updated) =
                match workspace.sub_detectors.entry(workspace_key) {
                    Entry::Occupied(entry) => {
                        let sub_detector = entry.into_mut();

                        if sub_detector.timestamp != new_timestamp {
                            let manifold = old_manifolds[sub_detector.manifold_id].take();
                            *sub_detector = SubDetector {
                                manifold_id: manifolds.len(),
                                timestamp: new_timestamp,
                                selected_contacts: 0,
                            };

                            manifolds.push(manifold);
                            (sub_detector, false)
                        } else {
                            (sub_detector, true)
                        }
                    }
                    Entry::Vacant(entry) => {
                        let sub_detector = SubDetector {
                            manifold_id: manifolds.len(),
                            selected_contacts: 0,
                            timestamp: new_timestamp,
                        };

                        let (id1, id2) = if flipped { (0, vid) } else { (vid, 0) };
                        manifolds.push(ContactManifold::with_data(
                            id1,
                            id2,
                            ManifoldData::default(),
                        ));

                        (entry.insert(sub_detector), false)
                    }
                };

            /*
             * Update the contact manifold if needed.
             */
            let manifold = &mut manifolds[sub_detector.manifold_id];

            if !manifold_updated {
                let mut canonical_mins1 = voxels1.voxel_center(key_low);
                let mut canonical_maxs1 = voxels1.voxel_center(key_high);

                for k in 0..DIM {
                    if key_low[k] != key_voxel[k] {
                        canonical_mins1[k] = canonical_mins1[k].max(domain2_1.mins[k]);
                    }

                    if key_high[k] != key_voxel[k] {
                        canonical_maxs1[k] = canonical_maxs1[k].min(domain2_1.maxs[k]);
                    }
                }

                let canonical_half_extents1 = (canonical_maxs1 - canonical_mins1) / 2.0 + radius1;
                let canonical_center1 = na::center(&canonical_mins1, &canonical_maxs1);

                let canonical_pseudo_cube1 = Cuboid::new(canonical_half_extents1);
                let canonical_pseudo_ball1 = RoundShape {
                    inner_shape: Cuboid::new(canonical_half_extents1 - radius1),
                    border_radius: radius1.x,
                };

                let canonical_shape1 = match voxels1.primitive_geometry() {
                    VoxelPrimitiveGeometry::PseudoBall { .. } => {
                        &canonical_pseudo_ball1 as &dyn Shape
                    }
                    VoxelPrimitiveGeometry::PseudoCube { .. } => {
                        &canonical_pseudo_cube1 as &dyn Shape
                    }
                };

                let canonical_pos12 = Translation::from(-canonical_center1) * pos12;

                // If we already computed contacts in the previous simulation step, their
                // local points are relative to the previously calculated canonical shape
                // which might not have the same local center as the one computed in this
                // step (because it’s based on the position of shape2 relative to voxels1).
                // So we need to adjust the local points to account for the position difference
                // and keep the point at the same "canonica-shape-space" location as in the previous frame.
                let prev_center = if flipped {
                    manifold
                        .subshape_pos2
                        .as_ref()
                        .map(|p| p.translation.vector)
                        .unwrap_or_default()
                } else {
                    manifold
                        .subshape_pos1
                        .as_ref()
                        .map(|p| p.translation.vector)
                        .unwrap_or_default()
                };
                let delta_center = canonical_center1.coords - prev_center;

                if flipped {
                    for pt in &mut manifold.points {
                        pt.local_p2 -= delta_center;
                    }
                } else {
                    for pt in &mut manifold.points {
                        pt.local_p1 -= delta_center;
                    }
                }

                // Update contacts.
                if flipped {
                    manifold.subshape_pos2 = Some(Isometry::from(canonical_center1));
                    let _ = dispatcher.contact_manifold_convex_convex(
                        &canonical_pos12.inverse(),
                        shape2,
                        canonical_shape1,
                        None,
                        None,
                        prediction,
                        manifold,
                    );
                } else {
                    manifold.subshape_pos1 = Some(Isometry::from(canonical_center1));
                    let _ = dispatcher.contact_manifold_convex_convex(
                        &canonical_pos12,
                        canonical_shape1,
                        shape2,
                        None,
                        None,
                        prediction,
                        manifold,
                    );
                }
            }

            /*
             * Filter-out points that don’t belong to this block.
             */
            let test_voxel = Cuboid::new(radius1 + Vector::repeat(1.0e-2));
            let penetration_dir1 = if flipped {
                manifold.local_n2
            } else {
                manifold.local_n1
            };

            for (i, pt) in manifold.points.iter().enumerate() {
                if pt.dist < 0.0 {
                    // If this is a penetration, double-check that we are not hitting the
                    // interior of the infinitely expanded canonical shape by checking if
                    // the opposite normal would have led to a better vector.
                    let cuboid1 = Cuboid::new(radius1);
                    let sp1 = cuboid1.local_support_point(&-penetration_dir1) + center1.coords;
                    let sm2 = shape2
                        .as_support_map()
                        .expect("Unsupported collision pair.");
                    let sp2 = sm2.support_point(pos12, &penetration_dir1);
                    let test_dist = (sp2 - sp1).dot(&-penetration_dir1);
                    let keep = test_dist < pt.dist;

                    if !keep {
                        // We don’t want to keep this point as it is an incorrect penetration
                        // caused by the canonical shape expansion.
                        continue;
                    }
                }

                let pt_in_voxel_space = if flipped {
                    manifold.subshape_pos2.transform_point(&pt.local_p2) - center1.coords
                } else {
                    manifold.subshape_pos1.transform_point(&pt.local_p1) - center1.coords
                };
                sub_detector.selected_contacts |=
                    (test_voxel.contains_local_point(&pt_in_voxel_space) as u32) << i;
            }
        }
    }

    // Remove contacts marked as ignored.
    for sub_detector in workspace.sub_detectors.values() {
        if sub_detector.timestamp == new_timestamp {
            let manifold = &mut manifolds[sub_detector.manifold_id];
            let mut k = 0;
            manifold.points.retain(|_| {
                let keep = (sub_detector.selected_contacts & (1 << k)) != 0;
                k += 1;
                keep
            });
        }
    }

    // Remove detectors which no longer index a valid manifold.
    workspace
        .sub_detectors
        .retain(|_, detector| detector.timestamp == new_timestamp);
}
