use crate::math::{Real, UnitVector, Vector, DEFAULT_EPSILON};
use crate::query::{CanonicalSplit, Split, SplitResult};
use crate::shape::{Segment, TriMesh};
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
            let idx = new_indices[k];
            for ic in 0..3 {
                let ia = (ic + 1) % 3;
                let ib = (ic + 1) % 3;
                let idx_a = idx[ia];
                let idx_b = idx[ib];
                let idx_c = idx[ic];

                if colors[idx_a as usize] + colors[idx_b as usize] == CROSSING_EDGE {
                    let intersection_idx = *intersections_found
                        .entry(SortedPair::new(idx_a, idx_b))
                        .or_insert_with(|| {
                            let segment = Segment::new(
                                new_vertices[idx_a as usize],
                                new_vertices[idx[ib] as usize],
                            );
                            // Intersect the segment with the plane.
                            if let SplitResult::Pair(left, right) =
                                segment.local_split(local_axis, bias, epsilon)
                            {
                                let new_pt = left.b;
                                new_vertices.push(new_pt);
                                colors.push(0);
                                new_vertices.len() - 1
                            } else {
                                unreachable!()
                            }
                        });

                    // Compute the indices of the two triangles.
                    let new_tri_a = [idx_a, intersection_idx as u32, idx_c];
                    let new_tri_b = [intersection_idx as u32, idx_b, idx_c];
                    // Replace the current triangle, and push the new one.
                    new_indices[k] = new_tri_a;
                    new_indices.push(new_tri_b);
                }
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

            if colors[idx[0]] == 1 || colors[idx[1]] == 1 || colors[idx[2]] == 1 {
                indices_lhs.push([remap[idx[0]].0, remap[idx[1]].0, remap[idx[2]].0]);
            } else if colors[idx[0]] == 2 || colors[idx[1]] == 2 || colors[idx[2]] == 2 {
                indices_rhs.push([remap[idx[0]].1, remap[idx[1]].1, remap[idx[2]].1]);
            } else {
                // The colors are all 0, so push into both trimeshes.
                indices_lhs.push([remap[idx[0]].0, remap[idx[1]].0, remap[idx[2]].0]);
                indices_rhs.push([remap[idx[0]].1, remap[idx[1]].1, remap[idx[2]].1]);
            }
        }

        let mesh_lhs = TriMesh::new(vertices_lhs, indices_lhs);
        let mesh_rhs = TriMesh::new(vertices_rhs, indices_rhs);
        SplitResult::Pair(mesh_lhs, mesh_rhs)
    }
}
