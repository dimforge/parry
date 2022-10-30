use crate::bounding_volume::{BoundingVolume, AABB};
use crate::math::{Isometry, Real};
use crate::query::contact_manifolds::contact_manifolds_workspace::{
    TypedWorkspaceData, WorkspaceData,
};
use crate::query::contact_manifolds::ContactManifoldsWorkspace;
use crate::query::query_dispatcher::PersistentQueryDispatcher;
use crate::query::ContactManifold;
use crate::shape::{Shape, TriMesh};
use crate::utils::hashmap::HashMap;

#[cfg_attr(feature = "serde-serialize", derive(Serialize, Deserialize))]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Deserialize, rkyv::Serialize)
)]
#[derive(Clone)]
pub struct TriMeshShapeContactManifoldsWorkspace {
    interferences: Vec<u32>,
    local_aabb2: AABB,
    old_interferences: Vec<u32>,
    delayed_manifolds: Vec<u32>,
    vertex_set: HashMap<u32, ()>,
}

impl TriMeshShapeContactManifoldsWorkspace {
    pub fn new() -> Self {
        Self {
            interferences: Vec::new(),
            local_aabb2: AABB::new_invalid(),
            old_interferences: Vec::new(),
            delayed_manifolds: Vec::new(),
            vertex_set: HashMap::default(),
        }
    }
}

/// Computes the contact manifold between a triangle-mesh an a shape, both represented as `Shape` trait-objects.
pub fn contact_manifolds_trimesh_shape_shapes<ManifoldData, ContactData>(
    dispatcher: &dyn PersistentQueryDispatcher<ManifoldData, ContactData>,
    pos12: &Isometry<Real>,
    shape1: &dyn Shape,
    shape2: &dyn Shape,
    prediction: Real,
    manifolds: &mut Vec<ContactManifold<ManifoldData, ContactData>>,
    workspace: &mut Option<ContactManifoldsWorkspace>,
) where
    ManifoldData: Default,
    ContactData: Default + Copy,
{
    if let Some(trimesh1) = shape1.as_trimesh() {
        contact_manifolds_trimesh_shape(
            dispatcher, pos12, trimesh1, shape2, prediction, manifolds, workspace, false,
        )
    } else if let Some(trimesh2) = shape2.as_trimesh() {
        contact_manifolds_trimesh_shape(
            dispatcher,
            &pos12.inverse(),
            trimesh2,
            shape1,
            prediction,
            manifolds,
            workspace,
            true,
        )
    }
}

fn ensure_workspace_exists(workspace: &mut Option<ContactManifoldsWorkspace>) {
    if workspace
        .as_mut()
        .and_then(|w| w.0.downcast_mut::<TriMeshShapeContactManifoldsWorkspace>())
        .is_some()
    {
        return;
    }

    *workspace = Some(ContactManifoldsWorkspace(Box::new(
        TriMeshShapeContactManifoldsWorkspace::new(),
    )));
}

/// Computes the contact manifold between a triangle-mesh and a shape.
pub fn contact_manifolds_trimesh_shape<ManifoldData, ContactData>(
    dispatcher: &dyn PersistentQueryDispatcher<ManifoldData, ContactData>,
    pos12: &Isometry<Real>,
    trimesh1: &TriMesh,
    shape2: &dyn Shape,
    prediction: Real,
    manifolds: &mut Vec<ContactManifold<ManifoldData, ContactData>>,
    workspace: &mut Option<ContactManifoldsWorkspace>,
    flipped: bool,
) where
    ManifoldData: Default,
    ContactData: Default + Copy,
{
    ensure_workspace_exists(workspace);
    let workspace: &mut TriMeshShapeContactManifoldsWorkspace =
        workspace.as_mut().unwrap().0.downcast_mut().unwrap();

    /*
     * Compute interferences.
     */
    // TODO: somehow precompute the AABB and reuse it?
    let mut new_local_aabb2 = shape2.compute_aabb(&pos12).loosened(prediction);
    let same_local_aabb2 = workspace.local_aabb2.contains(&new_local_aabb2);
    let mut old_manifolds = Vec::new();

    if !same_local_aabb2 {
        let extra_margin =
            (new_local_aabb2.maxs - new_local_aabb2.mins).map(|e| (e / 10.0).min(0.1));
        new_local_aabb2.mins -= extra_margin;
        new_local_aabb2.maxs += extra_margin;

        let local_aabb2 = new_local_aabb2; // .loosened(prediction * 2.0); // FIXME: what would be the best value?
        std::mem::swap(
            &mut workspace.old_interferences,
            &mut workspace.interferences,
        );

        std::mem::swap(manifolds, &mut old_manifolds);

        // This assertion may fire due to the invalid triangle_ids that the
        // near-phase may return (due to SIMD sentinels).
        //
        // assert_eq!(
        //     workspace
        //         .old_interferences
        //         .len()
        //         .min(trimesh1.num_triangles()),
        //     workspace.old_manifolds.len()
        // );

        workspace.interferences.clear();
        trimesh1
            .qbvh()
            .intersect_aabb(&local_aabb2, &mut workspace.interferences);
        workspace.local_aabb2 = local_aabb2;
    }

    /*
     * Dispatch to the specific solver by keeping the previous manifold if we already had one.
     */
    let new_interferences = &workspace.interferences;
    let mut old_inter_it = workspace.old_interferences.drain(..).peekable();
    let mut old_manifolds_it = old_manifolds.drain(..);

    // TODO: don't redispatch at each frame (we should probably do the same as
    // the heightfield).
    for (i, triangle_id) in new_interferences.iter().enumerate() {
        if *triangle_id >= trimesh1.num_triangles() as u32 {
            // Because of SIMD padding, the broad-phase may return triangle indices greater
            // than the max.
            continue;
        }

        if !same_local_aabb2 {
            loop {
                match old_inter_it.peek() {
                    Some(old_triangle_id) if *old_triangle_id < *triangle_id => {
                        let _ = old_inter_it.next();
                        let _ = old_manifolds_it.next();
                    }
                    _ => break,
                }
            }

            let manifold = if old_inter_it.peek() != Some(triangle_id) {
                let (id1, id2) = if flipped {
                    (0, *triangle_id)
                } else {
                    (*triangle_id, 0)
                };
                ContactManifold::with_data(id1, id2, ManifoldData::default())
            } else {
                // We already have a manifold for this triangle.
                let _ = old_inter_it.next();
                old_manifolds_it.next().unwrap()
            };

            manifolds.push(manifold);
        }

        let manifold = &mut manifolds[i];
        let triangle1 = trimesh1.triangle(*triangle_id);

        if flipped {
            let _ = dispatcher.contact_manifold_convex_convex(
                &pos12.inverse(),
                shape2,
                &triangle1,
                prediction,
                manifold,
            );
        } else {
            let _ = dispatcher
                .contact_manifold_convex_convex(pos12, &triangle1, shape2, prediction, manifold);
        }
    }

    /*
     *
     * Deal with internal edges.
     *
     */
    #[cfg(feature = "dim3")]
    {
        use crate::shape::FeatureId;
        use std::cmp::Ordering;

        // 1. Ingest all the face contacts.
        for (mid, manifold) in manifolds.iter().enumerate() {
            if let Some(deepest) = manifold.find_deepest_contact() {
                let tri_fid = if flipped { deepest.fid2 } else { deepest.fid1 };

                if tri_fid.is_face() {
                    let (tri_idx, tri) = if flipped {
                        (
                            trimesh1.indices()[manifold.subshape2 as usize],
                            trimesh1.triangle(manifold.subshape2),
                        )
                    } else {
                        (
                            trimesh1.indices()[manifold.subshape1 as usize],
                            trimesh1.triangle(manifold.subshape1),
                        )
                    };

                    let normal = if flipped {
                        manifold.local_n2
                    } else {
                        manifold.local_n1
                    };

                    if let Some(tri_normal) = tri.normal() {
                        // We check normal collinearity with an epsilon because sometimes,
                        // because of rounding errors, a contact may be identified as a face
                        // contact where it’s really just an edge contact.
                        if normal.dot(&tri_normal).abs() > 1.0 - 1.0e-4 {
                            let _ = workspace.vertex_set.insert(tri_idx[0], ());
                            let _ = workspace.vertex_set.insert(tri_idx[1], ());
                            let _ = workspace.vertex_set.insert(tri_idx[2], ());
                            // This as an actual face contact, continue without pushing to delayed_manifolds.
                            continue;
                        }
                    }
                }

                // NOTE: if we reach this line, then the contact isn’t an actual face contact.
                workspace.delayed_manifolds.push(mid as u32);
            }
        }

        // 2. Order by distance.
        workspace.delayed_manifolds.sort_by(|a, b| {
            let a = &manifolds[*a as usize];
            let b = &manifolds[*b as usize];

            let dist_a = a
                .find_deepest_contact()
                .map(|c| c.dist)
                .unwrap_or(Real::MAX);
            let dist_b = b
                .find_deepest_contact()
                .map(|c| c.dist)
                .unwrap_or(Real::MAX);

            dist_a.partial_cmp(&dist_b).unwrap_or(Ordering::Equal)
        });

        // 3. Deal with the edge/vertex contacts.
        for mid in &workspace.delayed_manifolds {
            let manifold = &mut manifolds[*mid as usize];

            let tri_idx = if flipped {
                trimesh1.indices()[manifold.subshape2 as usize]
            } else {
                trimesh1.indices()[manifold.subshape1 as usize]
            };

            manifold.points.retain(|pt| {
                let tri_fid = if flipped { pt.fid2 } else { pt.fid1 };

                match tri_fid.unpack() {
                    FeatureId::Face(_) => {
                        !workspace.vertex_set.contains_key(&tri_idx[0])
                            && !workspace.vertex_set.contains_key(&tri_idx[1])
                            && !workspace.vertex_set.contains_key(&tri_idx[2])
                    }
                    FeatureId::Edge(id) => {
                        !workspace.vertex_set.contains_key(&tri_idx[id as usize])
                            || !workspace
                                .vertex_set
                                .contains_key(&tri_idx[(id as usize + 1) % 3])
                    }
                    FeatureId::Vertex(id) => {
                        !workspace.vertex_set.contains_key(&tri_idx[id as usize])
                    }
                    FeatureId::Unknown => {
                        !workspace.vertex_set.contains_key(&tri_idx[0])
                            && !workspace.vertex_set.contains_key(&tri_idx[1])
                            && !workspace.vertex_set.contains_key(&tri_idx[2])
                    }
                }
            });

            // if !manifold.points.is_empty() {
            //     let normal = if flipped {
            //         manifold.local_n2
            //     } else {
            //         manifold.local_n1
            //     };
            //     println!("Keeping other contact: {:?}, {:?}", tri_idx, normal);
            // }

            let _ = workspace.vertex_set.insert(tri_idx[0], ());
            let _ = workspace.vertex_set.insert(tri_idx[1], ());
            let _ = workspace.vertex_set.insert(tri_idx[2], ());
        }

        workspace.vertex_set.clear();
        workspace.delayed_manifolds.clear();
    }
}

impl WorkspaceData for TriMeshShapeContactManifoldsWorkspace {
    fn as_typed_workspace_data(&self) -> TypedWorkspaceData {
        TypedWorkspaceData::TriMeshShapeContactManifoldsWorkspace(self)
    }

    fn clone_dyn(&self) -> Box<dyn WorkspaceData> {
        Box::new(self.clone())
    }
}
