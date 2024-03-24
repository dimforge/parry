use crate::query::ContactManifold;
use crate::shape::Triangle;
use crate::utils::hashmap::HashMap;

#[cfg(feature = "dim3")]
use crate::math::Real;

#[cfg_attr(feature = "serde-serialize", derive(Serialize, Deserialize))]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Deserialize, rkyv::Serialize),
    archive(check_bytes)
)]
#[derive(Default, Clone)]
#[allow(dead_code)] // We will need these for 2D too in the future.
pub struct InternalEdgesFixer {
    delayed_manifolds: Vec<u32>,
    vertex_set: HashMap<u32, ()>,
}

impl InternalEdgesFixer {
    #[cfg(feature = "dim2")]
    pub fn remove_invalid_contacts<ManifoldData, ContactData>(
        &mut self,
        _manifolds: &mut [ContactManifold<ManifoldData, ContactData>],
        _flipped: bool,
        _get_triangle: impl Fn(u32) -> Triangle,
        _get_triangle_indices: impl Fn(u32) -> [u32; 3],
    ) where
        ManifoldData: Default,
        ContactData: Default + Copy,
    {
        // Not yet implemented.
    }

    #[cfg(feature = "dim3")]
    pub fn remove_invalid_contacts<ManifoldData, ContactData>(
        &mut self,
        manifolds: &mut [ContactManifold<ManifoldData, ContactData>],
        flipped: bool,
        get_triangle: impl Fn(u32) -> Triangle,
        get_triangle_indices: impl Fn(u32) -> [u32; 3],
    ) where
        ManifoldData: Default,
        ContactData: Default + Copy,
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
                            get_triangle_indices(manifold.subshape2),
                            get_triangle(manifold.subshape2),
                        )
                    } else {
                        (
                            get_triangle_indices(manifold.subshape1),
                            get_triangle(manifold.subshape1),
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
                            let _ = self.vertex_set.insert(tri_idx[0], ());
                            let _ = self.vertex_set.insert(tri_idx[1], ());
                            let _ = self.vertex_set.insert(tri_idx[2], ());
                            // This as an actual face contact, continue without pushing to delayed_manifolds.
                            continue;
                        }
                    }
                }

                // NOTE: if we reach this line, then the contact isn’t an actual face contact.
                self.delayed_manifolds.push(mid as u32);
            }
        }

        // 2. Order by distance.
        self.delayed_manifolds.sort_by(|a, b| {
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
        for mid in &self.delayed_manifolds {
            let manifold = &mut manifolds[*mid as usize];

            let tri_idx = if flipped {
                get_triangle_indices(manifold.subshape2)
            } else {
                get_triangle_indices(manifold.subshape1)
            };

            manifold.points.retain(|pt| {
                let tri_fid = if flipped { pt.fid2 } else { pt.fid1 };

                match tri_fid.unpack() {
                    FeatureId::Face(_) => {
                        // This is a "false face" contact (a contact identified as face because
                        // of rounding errors but it’s probably just an edge contact.
                        // We don’t know exactly which edge it is, so we just ignore it
                        // if any of the vertices already exist in the set.
                        !self.vertex_set.contains_key(&tri_idx[0])
                            && !self.vertex_set.contains_key(&tri_idx[1])
                            && !self.vertex_set.contains_key(&tri_idx[2])
                    }
                    FeatureId::Edge(id) => {
                        !self.vertex_set.contains_key(&tri_idx[id as usize])
                            || !self
                                .vertex_set
                                .contains_key(&tri_idx[(id as usize + 1) % 3])
                    }
                    FeatureId::Vertex(id) => !self.vertex_set.contains_key(&tri_idx[id as usize]),
                    FeatureId::Unknown => {
                        // We don’t know the feature type here, so we just ignore it
                        // if any of the vertices already exist in the set.
                        !self.vertex_set.contains_key(&tri_idx[0])
                            && !self.vertex_set.contains_key(&tri_idx[1])
                            && !self.vertex_set.contains_key(&tri_idx[2])
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

            let _ = self.vertex_set.insert(tri_idx[0], ());
            let _ = self.vertex_set.insert(tri_idx[1], ());
            let _ = self.vertex_set.insert(tri_idx[2], ());
        }

        self.vertex_set.clear();
        self.delayed_manifolds.clear();
    }
}
