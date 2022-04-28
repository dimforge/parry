use crate::math::{Isometry, Point, Real, UnitVector, Vector};
use crate::query::visitors::BoundingVolumeIntersectionsVisitor;
use crate::query::{CanonicalSplit, PointQuery, Split, SplitResult};
use crate::shape::{Cuboid, FeatureId, Segment, Shape, TriMesh, Triangle};
use crate::utils::{hashmap::HashMap, SortedPair, WBasis};
use spade::{handles::FixedVertexHandle, ConstrainedDelaunayTriangulation, Triangulation as _};

impl CanonicalSplit for TriMesh {
    fn canonical_split(&self, axis: usize, bias: Real, epsilon: Real) -> SplitResult<Self> {
        // TODO optimize this.
        self.local_split(&Vector::ith_axis(axis), bias, epsilon)
    }
}

struct Triangulation {
    delaunay: ConstrainedDelaunayTriangulation<spade::Point2<Real>>,
    basis: [Vector<Real>; 2],
    basis_origin: Point<Real>,
    spade2index: HashMap<FixedVertexHandle, u32>,
    index2spade: HashMap<u32, FixedVertexHandle>,
}

impl Triangulation {
    fn new(axis: UnitVector<Real>, basis_origin: Point<Real>) -> Self {
        Triangulation {
            delaunay: ConstrainedDelaunayTriangulation::new(),
            basis: axis.orthonormal_basis(),
            basis_origin,
            spade2index: HashMap::default(),
            index2spade: HashMap::default(),
        }
    }

    fn project(&self, pt: Point<Real>) -> spade::Point2<Real> {
        let dpt = pt - self.basis_origin;
        spade::Point2::new(dpt.dot(&self.basis[0]), dpt.dot(&self.basis[1]))
    }

    fn add_edge(&mut self, id1: u32, id2: u32, points: &[Point<Real>]) {
        let proj1 = self.project(points[id1 as usize]);
        let proj2 = self.project(points[id2 as usize]);

        let handle1 = *self.index2spade.entry(id1).or_insert_with(|| {
            let h = self.delaunay.insert(proj1).unwrap();
            let _ = self.spade2index.insert(h, id1);
            h
        });

        let handle2 = *self.index2spade.entry(id2).or_insert_with(|| {
            let h = self.delaunay.insert(proj2).unwrap();
            let _ = self.spade2index.insert(h, id2);
            h
        });

        // NOTE: the naming of the `ConstrainedDelaunayTriangulation::can_add_constraint` method is misleading.
        if !self.delaunay.can_add_constraint(handle1, handle2) {
            let _ = self.delaunay.add_constraint(handle1, handle2);
        }
    }
}

impl Split for TriMesh {
    fn local_split(
        &self,
        local_axis: &UnitVector<Real>,
        bias: Real,
        epsilon: Real,
    ) -> SplitResult<Self> {
        let mut triangulation = if self.pseudo_normals().is_some() {
            Some(Triangulation::new(*local_axis, self.vertices()[0]))
        } else {
            None
        };

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
        let mut intersections_found = HashMap::default();
        let mut new_indices = indices.to_vec();
        let mut new_vertices = vertices.to_vec();

        for (tri_id, idx) in indices.iter().enumerate() {
            let mut idx = *idx;
            let mut intersection_features = (FeatureId::Unknown, FeatureId::Unknown);

            // First, find where the plane intersects the triangle.
            for ia in 0..3 {
                let ib = (ia + 1) % 3;
                let idx_a = idx[ia as usize];
                let idx_b = idx[ib as usize];

                let fid = match (colors[idx_a as usize], colors[idx_b as usize]) {
                    (1, 2) | (2, 1) => FeatureId::Edge(ia),
                    // NOTE: the case (_, 0) will be dealt with in the next loop iteration.
                    (0, _) => FeatureId::Vertex(ia),
                    _ => continue,
                };

                if intersection_features.0 == FeatureId::Unknown {
                    intersection_features.0 = fid;
                } else {
                    assert_eq!(intersection_features.1, FeatureId::Unknown);
                    intersection_features.1 = fid;
                }
            }

            // Helper that intersects an edge with the plane.
            let mut intersect_edge = |idx_a, idx_b| {
                *intersections_found
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
                            (new_vertices.len() - 1) as u32
                        } else {
                            unreachable!()
                        }
                    })
            };

            // Perform the intersection, push new triangles, and update
            // triangulation constraints if needed.
            match intersection_features {
                (_, FeatureId::Unknown) => {
                    // The plane doesn’t intersect the triangle, or intersects it at
                    // a single vertex, so we don’t have anything to do.
                    assert!(
                        matches!(intersection_features.0, FeatureId::Unknown)
                            || matches!(intersection_features.0, FeatureId::Vertex(_))
                    );
                }
                (FeatureId::Vertex(v1), FeatureId::Vertex(v2)) => {
                    // The plane intersects the triangle along one of its edge.
                    // We don’t have to split the triangle, but we need to add
                    // a constraint to the triangulation.
                    if let Some(triangulation) = &mut triangulation {
                        let id1 = idx[v1 as usize];
                        let id2 = idx[v2 as usize];
                        triangulation.add_edge(id1, id2, &new_vertices);
                    }
                }
                (FeatureId::Vertex(iv), FeatureId::Edge(ie))
                | (FeatureId::Edge(ie), FeatureId::Vertex(iv)) => {
                    // The plane splits the triangle into exactly two triangles.
                    let ia = ie;
                    let ib = (ie + 1) % 3;
                    let ic = (ie + 2) % 3;
                    let idx_a = idx[ia as usize];
                    let idx_b = idx[ib as usize];
                    let idx_c = idx[ic as usize];
                    assert_eq!(iv, ic);

                    let intersection_idx = intersect_edge(idx_a, idx_b);

                    // Compute the indices of the two triangles.
                    let new_tri_a = [idx_c, idx_a, intersection_idx as u32];
                    let new_tri_b = [idx_b, idx_c, intersection_idx as u32];

                    new_indices[tri_id] = new_tri_a;
                    new_indices.push(new_tri_b);

                    if let Some(triangulation) = &mut triangulation {
                        triangulation.add_edge(intersection_idx, idx_c, &new_vertices);
                    }
                }
                (FeatureId::Edge(mut e1), FeatureId::Edge(mut e2)) => {
                    // The plane splits the triangle into 1 + 2 triangles.
                    // First, make sure the edge indices are consecutive.
                    if e2 != (e1 + 1) % 3 {
                        std::mem::swap(&mut e1, &mut e2);
                    }

                    let ia = e2; // The first point of the second edge is the vertex shared by both edges.
                    let ib = (e2 + 1) % 3;
                    let ic = (e2 + 2) % 3;
                    let idx_a = idx[ia as usize];
                    let idx_b = idx[ib as usize];
                    let idx_c = idx[ic as usize];

                    let intersection1 = intersect_edge(idx_c, idx_a);
                    let intersection2 = intersect_edge(idx_a, idx_b);

                    let new_tri1 = [idx_a, intersection2, intersection1];
                    let new_tri2 = [intersection2, idx_b, idx_c];
                    let new_tri3 = [intersection2, idx_c, intersection1];
                    new_indices[tri_id] = new_tri1;
                    new_indices.push(new_tri2);
                    new_indices.push(new_tri3);

                    if let Some(triangulation) = &mut triangulation {
                        triangulation.add_edge(intersection1, intersection2, &new_vertices);
                    }
                }
                _ => unreachable!(),
            }
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

        // Push the triangulation if there is one.
        if let Some(triangulation) = triangulation {
            for face in triangulation.delaunay.inner_faces() {
                let vtx = face.vertices();
                let mut idx1 = [0; 3];
                let mut idx2 = [0; 3];
                for k in 0..3 {
                    let vid = triangulation.spade2index[&vtx[k].fix()];
                    assert_eq!(colors[vid as usize], 0);
                    idx1[k] = remap[vid as usize].0;
                    idx2[k] = remap[vid as usize].1;
                }

                let tri = Triangle::new(
                    vertices_lhs[idx1[0] as usize],
                    vertices_lhs[idx1[1] as usize],
                    vertices_lhs[idx1[2] as usize],
                );

                if self.contains_local_point(&tri.center()) {
                    indices_lhs.push(idx1);

                    idx2.swap(1, 2); // Flip orientation for the second half of the split.
                    indices_rhs.push(idx2);
                }
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
