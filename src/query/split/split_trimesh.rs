use crate::math::{Isometry, Real, UnitVector, Vector};
use crate::query::visitors::BoundingVolumeIntersectionsVisitor;
use crate::query::{CanonicalSplit, Split, SplitResult};
use crate::shape::{Cuboid, Segment, Shape, TriMesh};
use crate::utils::{hashmap::HashMap, SortedPair};

impl CanonicalSplit for TriMesh {
    fn canonical_split(&self, axis: usize, bias: Real, epsilon: Real) -> SplitResult<Self> {
        // TODO optimize this.
        self.local_split(&Vector::ith_axis(axis), bias, epsilon)
    }
}

impl Split for TriMesh {
    fn local_split(
        &self,
        local_axis: &UnitVector<Real>,
        bias: Real,
        epsilon: Real,
    ) -> SplitResult<Self> {
        // 1. Partition the vertices.
        let vertices = self.vertices();
        let indices = self.indices();
        let mut colors = vec![0u8; self.vertices().len()];

        // Color 0 = on plane.
        //       1 = on negative half-space.
        //       2 = on positive half-space.
        let mut found_negative = false;
        let mut found_positive = false;
        for (i, pt) in vertices.iter().enumerate() {
            let dist_to_plane = pt.coords.dot(&local_axis) - bias;
            if dist_to_plane < -epsilon {
                found_negative = true;
                colors[i] = 1;
            } else if dist_to_plane > epsilon {
                found_positive = true;
                colors[i] = 2;
            }
        }

        // Exit early if `self` isn’t crossed by the plane.
        if !found_negative {
            return SplitResult::Positive;
        }

        if !found_positive {
            return SplitResult::Negative;
        }

        // 2. Split the triangles.
        const CROSSING_EDGE: u8 = 3;
        let mut intersections_found = HashMap::default();
        let mut new_indices = indices.to_vec();
        let mut new_vertices = vertices.to_vec();
        let mut k = 0;

        while k != new_indices.len() {
            let mut ic = 0;
            while ic < 3 {
                let idx = new_indices[k];
                let ia = (ic + 1) % 3;
                let ib = (ic + 2) % 3;
                let idx_a = idx[ia];
                let idx_b = idx[ib];
                let idx_c = idx[ic];

                if colors[idx_a as usize] + colors[idx_b as usize] == CROSSING_EDGE {
                    let intersection_idx = *intersections_found
                        .entry(SortedPair::new(idx_a, idx_b))
                        .or_insert_with(|| {
                            let segment = Segment::new(
                                new_vertices[idx_a as usize],
                                new_vertices[idx_b as usize],
                            );
                            // Intersect the segment with the plane.
                            if let Some((intersection, _)) = segment
                                .local_split_and_get_intersection(local_axis, bias, epsilon)
                                .1
                            {
                                new_vertices.push(intersection);
                                colors.push(0);
                                new_vertices.len() - 1
                            } else {
                                unreachable!()
                            }
                        });

                    // Compute the indices of the two triangles.
                    let new_tri_a = [idx_c, idx_a, intersection_idx as u32];
                    let new_tri_b = [idx_b, idx_c, intersection_idx as u32];
                    // Replace the current triangle, and push the new one.
                    new_indices[k] = new_tri_a;
                    new_indices.push(new_tri_b);
                    // NOTE: we arranged the new triangle’s vertices such that, if
                    //       there is another intersection with `new_indices[k]`,
                    //       then that intersection can only happen with `ic == 2`
                    //       because we already know that the point at `idx[2]` lies
                    //       on the cutting plane.
                    ic = 2;
                    continue;
                }

                ic += 1;
            }

            k += 1;
        }

        // 3. Partition the new triangles into two trimeshes.
        let mut vertices_lhs = vec![];
        let mut vertices_rhs = vec![];
        let mut indices_lhs = vec![];
        let mut indices_rhs = vec![];
        let mut remap = vec![];

        for i in 0..new_vertices.len() {
            match colors[i] {
                0 => {
                    remap.push((vertices_lhs.len() as u32, vertices_rhs.len() as u32));
                    vertices_lhs.push(new_vertices[i]);
                    vertices_rhs.push(new_vertices[i]);
                }
                1 => {
                    remap.push((vertices_lhs.len() as u32, u32::MAX));
                    vertices_lhs.push(new_vertices[i]);
                }
                2 => {
                    remap.push((u32::MAX, vertices_rhs.len() as u32));
                    vertices_rhs.push(new_vertices[i]);
                }
                _ => unreachable!(),
            }
        }

        for idx in new_indices {
            let idx = [idx[0] as usize, idx[1] as usize, idx[2] as usize]; // Convert to usize.
            let colors = [colors[idx[0]], colors[idx[1]], colors[idx[2]]];
            let remap = [remap[idx[0]], remap[idx[1]], remap[idx[2]]];

            if colors[0] == 1 || colors[1] == 1 || colors[2] == 1 {
                assert!(colors[0] != 2 && colors[1] != 2 && colors[2] != 2);
                indices_lhs.push([remap[0].0, remap[1].0, remap[2].0]);
            } else if colors[0] == 2 || colors[1] == 2 || colors[2] == 2 {
                assert!(colors[0] != 1 && colors[1] != 1 && colors[2] != 1);
                indices_rhs.push([remap[0].1, remap[1].1, remap[2].1]);
            } else {
                // The colors are all 0, so push into both trimeshes.
                indices_lhs.push([remap[0].0, remap[1].0, remap[2].0]);
                indices_rhs.push([remap[0].1, remap[1].1, remap[2].1]);
            }
        }

        let mesh_lhs = TriMesh::new(vertices_lhs, indices_lhs);
        let mesh_rhs = TriMesh::new(vertices_rhs, indices_rhs);
        SplitResult::Pair(mesh_lhs, mesh_rhs)
    }

    fn intersection_with_local_cuboid(
        &self,
        cuboid: &Cuboid,
        cuboid_position: &Isometry<Real>,
        epsilon: Real,
    ) -> Option<Self> {
        let cuboid_aabb = cuboid.compute_aabb(cuboid_position);
        let mut intersecting_tris = vec![];
        let mut visitor = BoundingVolumeIntersectionsVisitor::new(&cuboid_aabb, |id| {
            intersecting_tris.push(*id);
            true
        });
        self.qbvh().traverse_depth_first(&mut visitor);

        if intersecting_tris.is_empty() {
            return None;
        }

        // First, very naive version that outputs a triangle soup without
        // index buffer (shared vertices are duplicated).
        let vertices = self.vertices();
        let indices = self.indices();

        let mut clip_workspace = vec![];
        let mut new_vertices = vec![];
        let mut new_indices = vec![];
        let aabb = cuboid.local_aabb();
        let inv_pos = cuboid_position.inverse();
        let mut to_clip = vec![];

        for tri in intersecting_tris {
            let idx = indices[tri as usize];
            to_clip.extend_from_slice(&[
                inv_pos * vertices[idx[0] as usize],
                inv_pos * vertices[idx[1] as usize],
                inv_pos * vertices[idx[2] as usize],
            ]);

            // There is no need to clip if the triangle is fully inside of the AABB.
            // Note that we can’t take a shortcut for the case where all the vertices are
            // outside of the AABB, because the AABB can still instersect the edges or face.
            if !(aabb.contains_local_point(&to_clip[0])
                && aabb.contains_local_point(&to_clip[1])
                && aabb.contains_local_point(&to_clip[2]))
            {
                aabb.clip_polygon_with_workspace(&mut to_clip, &mut clip_workspace);
            }

            if to_clip.len() >= 3 {
                let base_i = new_vertices.len();
                for i in 1..to_clip.len() - 1 {
                    new_indices.push([base_i as u32, (base_i + i) as u32, (base_i + i + 1) as u32]);
                }
                new_vertices.append(&mut to_clip);
            }
        }

        // The clipping outputs points in the local-space of the cuboid.
        // So we need to transform it back.
        for pt in &mut new_vertices {
            *pt = cuboid_position * *pt;
        }

        if new_vertices.len() >= 3 {
            Some(TriMesh::new(new_vertices, new_indices))
        } else {
            None
        }
    }
}
