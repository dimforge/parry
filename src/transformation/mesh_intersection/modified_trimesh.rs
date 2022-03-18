use crate::math::{Isometry, Point, Real};
use crate::shape::trimesh::{TopoFace, TopoHalfEdge, TopoVertex, TriMeshTopology};
use crate::shape::{TriMesh, Triangle};
use std::collections::{HashMap, HashSet};

#[derive(Copy, Clone, Debug)]
pub struct IntersectionHalfEdge {
    pub first_vertex: u32,
    pub second_vertex: u32,
    pub eid: u32,
}

impl IntersectionHalfEdge {
    pub fn new() -> Self {
        Self {
            first_vertex: u32::MAX,
            second_vertex: u32::MAX,
            eid: u32::MAX,
        }
    }
}

pub struct ModifiedTriMesh<'a> {
    pub tag: u8, // NOTE: only for debugging.
    pub trimesh: &'a TriMesh,
    // NOTE: the added vertices are in the local-space of the first trimesh.
    pub added_vertices: Vec<Point<Real>>,
    pub added_faces: Vec<[u32; 3]>,
    // The original triangle this face is split from.
    pub original_faces: Vec<u32>,
    pub added_topo: TriMeshTopology,
    pub half_edge_substitutes: HashMap<u32, u32>,
    // TODO: we need a hash-map only for the first time a face
    //       is split. If it’s split more than once, the subsequent splits
    //       could be stored in a Vec with direct indexed access.
    pub face_splits: HashMap<u32, u32>,
    pub edge_splits: HashSet<u32>, // For debugging
}

impl<'a> ModifiedTriMesh<'a> {
    pub fn new(trimesh: &'a TriMesh, tag: u8) -> Self {
        Self {
            tag,
            trimesh,
            added_vertices: vec![],
            added_faces: vec![],
            original_faces: vec![],
            added_topo: TriMeshTopology::default(),
            half_edge_substitutes: HashMap::default(),
            face_splits: HashMap::default(),
            edge_splits: HashSet::default(),
        }
    }

    // NOTE: this is mostly used for debugging, in order to visualize the whay faces have been cut.
    pub fn into_trimesh(mut self, pos: &Isometry<Real>) -> TriMesh {
        let new_indices: Vec<_> = self
            .trimesh
            .indices()
            .into_iter()
            .copied()
            .chain(self.added_faces.into_iter())
            .enumerate()
            .filter(|(i, idx)| !self.face_splits.contains_key(&(*i as u32)) && idx[0] != u32::MAX)
            .map(|(_, idx)| idx)
            .collect();

        let mut new_vertices: Vec<_> = self.trimesh.vertices().iter().map(|pt| pos * pt).collect();
        new_vertices.append(&mut self.added_vertices);

        TriMesh::new(new_vertices, new_indices)
    }

    pub fn resolve_triangles(
        &self,
        id: u32,
        pos: &Isometry<Real>,
        out: &mut Vec<(u32, Triangle, [u32; 3])>,
    ) {
        if let Some(face_splits) = self.face_splits.get(&id) {
            self.resolve_triangles(*face_splits, pos, out);
            self.resolve_triangles(*face_splits + 1, pos, out);
            self.resolve_triangles(*face_splits + 2, pos, out);
        } else if let Some((tri, idx)) = self.face_vertices(id, pos) {
            out.push((id, tri, idx));
        }
    }

    pub fn vertex(&self, vid: u32, pos: &Isometry<Real>) -> Point<Real> {
        let vertices = self.trimesh.vertices();
        if (vid as usize) < vertices.len() {
            pos * vertices[vid as usize]
        } else {
            self.added_vertices[vid as usize - vertices.len()]
        }
    }

    pub fn original_face_id(&self, fid: u32) -> u32 {
        let base_len = self.trimesh.indices().len();
        if fid >= base_len as u32 {
            self.original_faces[fid as usize - base_len]
        } else {
            fid
        }
    }

    pub fn original_face_vertices(&self, fid: u32, pos: &Isometry<Real>) -> Triangle {
        self.trimesh
            .triangle(self.original_face_id(fid))
            .transformed(pos)
    }

    pub fn face_vertices(&self, fid: u32, pos: &Isometry<Real>) -> Option<(Triangle, [u32; 3])> {
        let idx = self.face_indices(fid);
        // If this is not an invalid triangle (added by edge-splitting), push it.
        if idx[0] != u32::MAX {
            Some((
                Triangle::new(
                    self.vertex(idx[0], pos),
                    self.vertex(idx[1], pos),
                    self.vertex(idx[2], pos),
                ),
                idx,
            ))
        } else {
            None
        }
    }

    pub fn half_edge_base_id(&self) -> u32 {
        (self.trimesh.topology().half_edges.len() + self.added_topo.half_edges.len()) as u32
    }

    pub fn vertex_base_id(&self) -> u32 {
        (self.trimesh.vertices().len() + self.added_vertices.len()) as u32
    }

    pub fn face_base_id(&self) -> u32 {
        (self.trimesh.topology().faces.len() + self.added_topo.faces.len()) as u32
    }

    pub fn face_half_edges_ids(&self, fid: u32) -> [u32; 3] {
        assert!(!self.face_splits.contains_key(&fid));
        let prev_face_len = self.trimesh.topology().faces.len();
        let first_half_edge = if let Some(face) = self.trimesh.topology().faces.get(fid as usize) {
            face.half_edge
        } else {
            self.added_topo.faces[fid as usize - prev_face_len].half_edge
        };

        let mut result = [first_half_edge; 3];
        for k in 1..3 {
            let half_edge = self.half_edge(result[k - 1]);
            result[k] = half_edge.next;
        }
        result
    }

    pub fn face_indices_and_edge_index(&self, half_edge: &TopoHalfEdge) -> ([u32; 3], u8) {
        let indices = self.face_indices(half_edge.face);
        for i in 0..3 {
            if indices[i] == half_edge.vertex {
                return (indices, i as u8);
            }
        }

        unreachable!();
    }

    pub fn face_indices(&self, fid: u32) -> [u32; 3] {
        assert!(!self.face_splits.contains_key(&fid));
        if let Some(idx) = self.trimesh.indices().get(fid as usize) {
            *idx
        } else {
            let prev_face_len = self.trimesh.topology().faces.len();
            self.added_faces[fid as usize - prev_face_len]
        }
    }

    pub fn half_edge(&self, id: u32) -> &TopoHalfEdge {
        let id = self.half_edge_substitutes.get(&id).copied().unwrap_or(id);
        if let Some(hedge) = self.trimesh.topology().half_edges.get(id as usize) {
            hedge
        } else {
            let prev_half_edges = self.trimesh.topology().half_edges.len();
            &self.added_topo.half_edges[id as usize - prev_half_edges]
        }
    }

    pub fn substitute_half_edge(&mut self, id_to_substitute: u32) -> (&mut TopoHalfEdge, u32) {
        let mut id = self
            .half_edge_substitutes
            .get(&id_to_substitute)
            .copied()
            .unwrap_or(id_to_substitute);
        let base_len = self.trimesh.topology().half_edges.len();

        if let Some(hedge) = self.trimesh.topology().half_edges.get(id as usize) {
            id = self.half_edge_base_id();
            self.added_topo.half_edges.push(*hedge);
            let _ = self.half_edge_substitutes.insert(id_to_substitute, id);
        }

        (&mut self.added_topo.half_edges[id as usize - base_len], id)
    }

    pub fn assert_face_isnt_degenerate(&self, idx: [u32; 3], pos: &Isometry<Real>) {
        let tri = Triangle::new(
            self.vertex(idx[0], pos),
            self.vertex(idx[1], pos),
            self.vertex(idx[2], pos),
        );
        if tri.area() == 0.0 {
            // dbg!(tri);
        }
        // assert_ne!(tri.area(), 0.0);
    }

    pub fn push_original_face(&mut self, fid: u32) {
        let base_len = self.trimesh.indices().len();
        if (fid as usize) < base_len {
            self.original_faces.push(fid);
        } else {
            self.original_faces
                .push(self.original_faces[fid as usize - base_len]);
        }
    }

    pub fn split_edge(
        &mut self,
        half_edge_id: u32,
        new_pt: Point<Real>,
        intersection: &mut IntersectionHalfEdge,
        pos: &Isometry<Real>,
    ) -> [u32; 2] {
        // println!("Splitting half-egde: {:?}", half_edge_id);
        let vertex_base_id = self.vertex_base_id();
        let face_base_id = self.face_base_id();
        let half_edge_base_id = self.half_edge_base_id();
        let new_pt_id = self.vertex_base_id();

        let half_edge = *self.half_edge(half_edge_id);
        let fid = half_edge.face;
        let (face_indices, edge) = self.face_indices_and_edge_index(&half_edge);

        let twin_half_edge = *self.half_edge(half_edge.twin);
        let twin_fid = twin_half_edge.face;
        let (twin_face_indices, twin_edge) = self.face_indices_and_edge_index(&twin_half_edge);

        assert!(self.edge_splits.insert(half_edge_id));
        assert!(self.edge_splits.insert(half_edge.twin));

        // Grab the indices of the half-edges that will survive the split (one half-edge per triangle after-split).
        let next_half_edge = half_edge.next;
        let next_next_half_edge = self.half_edge(next_half_edge).next;

        let twin_next_half_edge = twin_half_edge.next;
        let twin_next_next_half_edge = self.half_edge(twin_next_half_edge).next;

        // Add 1 vertex.
        self.added_vertices.push(new_pt);
        self.added_topo.vertices.push(TopoVertex {
            half_edge: half_edge_base_id,
        });

        // Add 4 faces.
        assert_eq!(
            face_indices[(edge as usize + 0) % 3],
            twin_face_indices[(twin_edge as usize + 1) % 3]
        );
        assert_eq!(
            face_indices[(edge as usize + 1) % 3],
            twin_face_indices[(twin_edge as usize + 0) % 3]
        );

        let new_face_indices = [
            [
                face_indices[(edge as usize + 2) % 3],
                face_indices[(edge as usize + 0) % 3],
                new_pt_id,
            ],
            [
                face_indices[(edge as usize + 1) % 3],
                face_indices[(edge as usize + 2) % 3],
                new_pt_id,
            ],
            [u32::MAX; 3], // We always split triangles into 3 triangles. So, for the third tri, use invalid indices.
        ];
        let new_twin_face_indices = [
            [
                twin_face_indices[(twin_edge as usize + 2) % 3],
                twin_face_indices[(twin_edge as usize + 0) % 3],
                new_pt_id,
            ],
            [
                twin_face_indices[(twin_edge as usize + 1) % 3],
                twin_face_indices[(twin_edge as usize + 2) % 3],
                new_pt_id,
            ],
            [u32::MAX; 3],
        ];

        // println!("[{}] Splitting edge face {fid}", self.tag);
        // println!("[{}] Splitting edge face {twin_fid}", self.tag);
        assert!(
            self.face_splits.insert(fid, face_base_id).is_none(),
            "Found duplicate face split."
        );
        assert!(
            self.face_splits
                .insert(twin_fid, face_base_id + 3)
                .is_none(),
            "Found duplicate face split."
        );

        // self.assert_face_isnt_degenerate(new_face_indices[0], pos);
        // self.assert_face_isnt_degenerate(new_face_indices[1], pos);
        // self.assert_face_isnt_degenerate(new_twin_face_indices[0], pos);
        // self.assert_face_isnt_degenerate(new_twin_face_indices[1], pos);

        for k in 0..3 {
            self.added_faces.push(new_face_indices[k]);
            self.push_original_face(fid);
        }
        for k in 0..3 {
            self.added_faces.push(new_twin_face_indices[k]);
            self.push_original_face(twin_fid);
        }
        self.added_topo.faces.push(TopoFace {
            half_edge: next_next_half_edge,
        });
        self.added_topo.faces.push(TopoFace {
            half_edge: next_half_edge,
        });
        self.added_topo.faces.push(TopoFace {
            half_edge: u32::MAX,
        });
        self.added_topo.faces.push(TopoFace {
            half_edge: twin_next_next_half_edge,
        });
        self.added_topo.faces.push(TopoFace {
            half_edge: twin_next_half_edge,
        });
        self.added_topo.faces.push(TopoFace {
            half_edge: u32::MAX,
        });

        // Add the 8 new half-edges (we add them face-after-face, in CCW order).
        self.added_topo.half_edges.push(TopoHalfEdge {
            next: half_edge_base_id + 1,
            twin: half_edge_base_id + 7,
            vertex: new_face_indices[0][1],
            face: face_base_id + 0,
        });
        self.added_topo.half_edges.push(TopoHalfEdge {
            next: next_next_half_edge,
            twin: half_edge_base_id + 2,
            vertex: new_pt_id,
            face: face_base_id + 0,
        });
        self.added_topo.half_edges.push(TopoHalfEdge {
            next: half_edge_base_id + 3,
            twin: half_edge_base_id + 1,
            vertex: new_face_indices[1][1],
            face: face_base_id + 1,
        });
        self.added_topo.half_edges.push(TopoHalfEdge {
            next: next_half_edge,
            twin: half_edge_base_id + 4,
            vertex: new_pt_id,
            face: face_base_id + 1,
        });
        self.added_topo.half_edges.push(TopoHalfEdge {
            next: half_edge_base_id + 5,
            twin: half_edge_base_id + 3,
            vertex: new_twin_face_indices[0][1],
            face: face_base_id + 3,
        });
        self.added_topo.half_edges.push(TopoHalfEdge {
            next: twin_next_next_half_edge,
            twin: half_edge_base_id + 6,
            vertex: new_pt_id,
            face: face_base_id + 3,
        });
        self.added_topo.half_edges.push(TopoHalfEdge {
            next: half_edge_base_id + 7,
            twin: half_edge_base_id + 5,
            vertex: new_twin_face_indices[1][1],
            face: face_base_id + 4,
        });
        self.added_topo.half_edges.push(TopoHalfEdge {
            next: twin_next_half_edge,
            twin: half_edge_base_id + 0,
            vertex: new_pt_id,
            face: face_base_id + 4,
        });

        // Modify the surviving adjascent half-edges:
        // - update `next`.
        // - update `face`.
        let (next_next, _) = self.substitute_half_edge(next_next_half_edge);
        next_next.next = half_edge_base_id + 0;
        next_next.face = face_base_id + 0;

        let (next, _) = self.substitute_half_edge(next_half_edge);
        next.next = half_edge_base_id + 2;
        next.face = face_base_id + 1;

        let (twin_next_next, _) = self.substitute_half_edge(twin_next_next_half_edge);
        twin_next_next.next = half_edge_base_id + 4;
        twin_next_next.face = face_base_id + 3;

        let (twin_next, _) = self.substitute_half_edge(twin_next_half_edge);
        twin_next.next = half_edge_base_id + 6;
        twin_next.face = face_base_id + 4;

        if intersection.first_vertex == u32::MAX {
            intersection.first_vertex = vertex_base_id;
        } else {
            intersection.second_vertex = vertex_base_id;
            let base = self.added_topo.half_edges.len() - 8;
            for k in 0..8 {
                // if self.added_topo.half_edges[base + k].vertex == intersection.first_vertex {
                if self.half_edge(half_edge_base_id + k).vertex == intersection.first_vertex {
                    intersection.eid = half_edge_base_id + k as u32;
                }
            }
        }

        [half_edge_base_id, half_edge_base_id + 3]
    }

    pub fn split_face(
        &mut self,
        fid: u32,
        new_pt: Point<Real>,
        intersection: &mut IntersectionHalfEdge,
        pos: &Isometry<Real>,
    ) -> ([u32; 3], [u32; 3]) {
        // println!("[{}] Splitting face {fid}", self.tag);
        let vertex_base_id = self.vertex_base_id();
        let face_base_id = self.face_base_id();
        let half_edge_base_id = self.half_edge_base_id();
        let idx = self.face_indices(fid);
        let new_pt_id = self.vertex_base_id();
        let half_edge_ids = self.face_half_edges_ids(fid);

        // Add 1 vertex.
        self.added_vertices.push(new_pt);
        self.added_topo.vertices.push(TopoVertex {
            half_edge: half_edge_base_id,
        });

        // Add 3 faces.
        assert!(
            self.face_splits.insert(fid, face_base_id).is_none(),
            "Found duplicate face split."
        );

        let tri1 = [idx[0], idx[1], new_pt_id];
        let tri2 = [idx[1], idx[2], new_pt_id];
        let tri3 = [idx[2], idx[0], new_pt_id];
        self.added_faces.push(tri1);
        self.added_faces.push(tri2);
        self.added_faces.push(tri3);
        self.push_original_face(fid);
        self.push_original_face(fid);
        self.push_original_face(fid);

        self.assert_face_isnt_degenerate(tri1, pos);
        self.assert_face_isnt_degenerate(tri2, pos);
        self.assert_face_isnt_degenerate(tri3, pos);

        self.added_topo.faces.push(TopoFace {
            half_edge: half_edge_ids[0],
        });
        self.added_topo.faces.push(TopoFace {
            half_edge: half_edge_ids[1],
        });
        self.added_topo.faces.push(TopoFace {
            half_edge: half_edge_ids[2],
        });

        // Add 3 half-edges starting at new_point.
        for k in 0u32..3 {
            self.added_topo.half_edges.push(TopoHalfEdge {
                next: half_edge_ids[k as usize],
                twin: half_edge_base_id + k + 3,
                vertex: new_pt_id,
                face: face_base_id + k,
            });
        }

        // Add 3 half-edges ending at new_point.
        for k in 0u32..3 {
            self.added_topo.half_edges.push(TopoHalfEdge {
                next: half_edge_base_id + (k + 2) % 3,
                twin: half_edge_base_id + k,
                vertex: idx[k as usize],
                face: face_base_id + (k + 2) % 3,
            });
        }

        // From existing half-edges:
        // - Modify the `next` index.
        // - Modify the `face` index.
        let nexts = [4, 5, 3];
        for k in 0u32..3 {
            let (hedge, _) = self.substitute_half_edge(half_edge_ids[k as usize]);
            hedge.face = face_base_id + k;
            hedge.next = half_edge_base_id + nexts[k as usize];
        }

        if intersection.first_vertex == u32::MAX {
            intersection.first_vertex = vertex_base_id;
        } else {
            intersection.second_vertex = vertex_base_id;
            // let base = self.added_topo.half_edges.len() - 6;
            for k in 0..6 {
                // if self.added_topo.half_edges[base + k].vertex == intersection.first_vertex {
                if self.half_edge(half_edge_base_id + k).vertex == intersection.first_vertex {
                    assert!(k >= 3);
                    intersection.eid = half_edge_base_id + k as u32;
                }
            }
        }

        (
            [face_base_id, face_base_id + 1, face_base_id + 2],
            [
                half_edge_base_id,
                half_edge_base_id + 1,
                half_edge_base_id + 2,
            ],
        )
    }

    pub fn insert_vertex(&mut self, triangle: u32, vertex: Point<Real>) -> u32 {
        todo!()
    }

    pub fn insert_edge(&mut self, triangle: u32, a: u32, b: u32) -> u32 {
        todo!()
    }
}
