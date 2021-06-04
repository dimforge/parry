use crate::math::{Point, Real, Vector, DIM};
use crate::shape::{FeatureId, PolygonalFeature, PolygonalFeatureMap, SupportMap};
// use crate::transformation;
use crate::utils::hashmap::{Entry, HashMap};
use crate::utils::{self, SortedPair};
use na::{self, ComplexField, Point2, Unit};
use std::f64;

#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[derive(PartialEq, Debug, Copy, Clone)]
pub struct Vertex {
    pub first_adj_face_or_edge: u32,
    pub num_adj_faces_or_edge: u32,
}

#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[derive(PartialEq, Debug, Copy, Clone)]
pub struct Edge {
    pub vertices: Point2<u32>,
    pub faces: Point2<u32>,
    pub dir: Unit<Vector<Real>>,
    deleted: bool,
}

impl Edge {
    fn other_triangle(&self, id: u32) -> u32 {
        if id == self.faces[0] {
            self.faces[1]
        } else {
            self.faces[0]
        }
    }
}

#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[derive(PartialEq, Debug, Copy, Clone)]
pub struct Face {
    pub first_vertex_or_edge: u32,
    pub num_vertices_or_edges: u32,
    pub normal: Unit<Vector<Real>>,
}

#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[derive(PartialEq, Debug, Copy, Clone)]
struct Triangle {
    vertices: [u32; 3],
    edges: [u32; 3],
    normal: Vector<Real>,
    parent_face: Option<u32>,
    is_degenerate: bool,
}

impl Triangle {
    fn next_edge_id(&self, id: u32) -> u32 {
        for i in 0..3 {
            if self.edges[i] == id {
                return (i as u32 + 1) % 3;
            }
        }

        unreachable!()
    }
}

#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[derive(PartialEq, Debug, Clone)]
/// A convex polyhedron without degenerate faces.
pub struct ConvexPolyhedron {
    points: Vec<Point<Real>>,
    vertices: Vec<Vertex>,
    faces: Vec<Face>,
    edges: Vec<Edge>,
    // Faces adjascent to a vertex.
    faces_adj_to_vertex: Vec<u32>,
    // Edges adjascent to a vertex.
    edges_adj_to_vertex: Vec<u32>,
    // Edges adjascent to a face.
    edges_adj_to_face: Vec<u32>,
    // Vertices adjascent to a face.
    vertices_adj_to_face: Vec<u32>,
}

impl ConvexPolyhedron {
    /// Creates a new convex polyhedron from an arbitrary set of points.
    ///
    /// This explicitly computes the convex hull of the given set of points. Use
    /// Returns `None` if the convex hull computation failed.
    pub fn from_convex_hull(points: &[Point<Real>]) -> Option<ConvexPolyhedron> {
        let (vertices, indices) = crate::transformation::convex_hull(points);
        Self::from_convex_mesh(vertices, &indices)
    }

    /// Attempts to create a new solid assumed to be convex from the set of points and indices.
    ///
    /// The given points and index information are assumed to describe a convex polyhedron.
    /// It it is not, weird results may be produced.
    ///
    /// # Return
    ///
    /// Retruns `None` if he given solid is not manifold (contains t-junctions, not closed, etc.)
    pub fn from_convex_mesh(
        points: Vec<Point<Real>>,
        indices: &[[u32; DIM]],
    ) -> Option<ConvexPolyhedron> {
        let eps = ComplexField::sqrt(crate::math::DEFAULT_EPSILON);

        let mut vertices = Vec::new();
        let mut edges = Vec::<Edge>::new();
        let mut faces = Vec::<Face>::new();
        let mut triangles = Vec::new();
        let mut edge_map = HashMap::default();

        let mut faces_adj_to_vertex = Vec::new();
        let mut edges_adj_to_vertex = Vec::new();
        let mut edges_adj_to_face = Vec::new();
        let mut vertices_adj_to_face = Vec::new();

        //// Euler characteristic.
        let nedges = points.len() + indices.len() - 2;
        edges.reserve(nedges);

        /*
         *  Initialize triangles and edges adjacency information.
         */
        for vtx in indices {
            let mut edges_id = [u32::MAX; DIM];
            let face_id = triangles.len();

            assert!(vtx[0] != vtx[1]);
            assert!(vtx[0] != vtx[2]);
            assert!(vtx[2] != vtx[1]);

            for i1 in 0..3 {
                // Deal with edges.
                let i2 = (i1 + 1) % 3;
                let key = SortedPair::new(vtx[i1], vtx[i2]);

                match edge_map.entry(key) {
                    Entry::Occupied(e) => {
                        let edge = &mut edges[*e.get() as usize];
                        let out_face_id = &mut edge.faces[1];

                        if *out_face_id == u32::MAX {
                            edges_id[i1] = *e.get();
                            *out_face_id = face_id as u32
                        } else {
                            // We have a t-junction.
                            return None;
                        }
                    }
                    Entry::Vacant(e) => {
                        edges_id[i1] = *e.insert(edges.len() as u32);

                        let dir = Unit::try_new(
                            points[vtx[i2] as usize] - points[vtx[i1] as usize],
                            crate::math::DEFAULT_EPSILON,
                        );

                        edges.push(Edge {
                            vertices: Point2::new(vtx[i1], vtx[i2]),
                            faces: Point2::new(face_id as u32, u32::MAX),
                            dir: dir.unwrap_or(Vector::x_axis()),
                            deleted: dir.is_none(),
                        });
                    }
                }
            }

            let normal = utils::ccw_face_normal([
                &points[vtx[0] as usize],
                &points[vtx[1] as usize],
                &points[vtx[2] as usize],
            ]);
            let triangle = Triangle {
                vertices: *vtx,
                edges: edges_id,
                normal: normal.map(|n| *n).unwrap_or(Vector::zeros()),
                parent_face: None,
                is_degenerate: normal.is_none(),
            };

            triangles.push(triangle);
        }

        // Find edges that must be deleted.

        for e in &mut edges {
            let tri1 = triangles.get(e.faces[0] as usize)?;
            let tri2 = triangles.get(e.faces[1] as usize)?;
            if tri1.normal.dot(&tri2.normal) > 1.0 - eps {
                e.deleted = true;
            }
        }

        /*
         * Extract faces by following  contours.
         */
        for i in 0..triangles.len() {
            if triangles[i].parent_face.is_none() {
                for j1 in 0..3 {
                    if !edges[triangles[i].edges[j1] as usize].deleted {
                        // Create a new face, setup its first edge/vertex and construct it.
                        let new_face_id = faces.len();
                        let mut new_face = Face {
                            first_vertex_or_edge: edges_adj_to_face.len() as u32,
                            num_vertices_or_edges: 1,
                            normal: Unit::new_unchecked(triangles[i].normal),
                        };

                        edges_adj_to_face.push(triangles[i].edges[j1]);
                        vertices_adj_to_face.push(triangles[i].vertices[j1]);

                        let j2 = (j1 + 1) % 3;
                        let start_vertex = triangles[i].vertices[j1];

                        // NOTE: variables ending with _id are identifier on the
                        // fields of a triangle. Other variables are identifier on
                        // the triangles/edges/vertices arrays.
                        let mut curr_triangle = i;
                        let mut curr_edge_id = j2;

                        while triangles[curr_triangle].vertices[curr_edge_id] != start_vertex {
                            let curr_edge = triangles[curr_triangle].edges[curr_edge_id];
                            let curr_vertex = triangles[curr_triangle].vertices[curr_edge_id];
                            // NOTE: we should use this assertion. However, it can currently
                            // happen if there are some isolated non-deleted edges due to
                            // rounding errors.
                            //
                            // assert!(triangles[curr_triangle].parent_face.is_none());
                            triangles[curr_triangle].parent_face = Some(new_face_id as u32);

                            if !edges[curr_edge as usize].deleted {
                                edges_adj_to_face.push(curr_edge);
                                vertices_adj_to_face.push(curr_vertex);
                                new_face.num_vertices_or_edges += 1;

                                curr_edge_id = (curr_edge_id + 1) % 3;
                            } else {
                                // Find adjacent edge on the next triangle.
                                curr_triangle = edges[curr_edge as usize]
                                    .other_triangle(curr_triangle as u32)
                                    as usize;
                                curr_edge_id =
                                    triangles[curr_triangle].next_edge_id(curr_edge) as usize;
                                assert!(
                                    triangles[curr_triangle].vertices[curr_edge_id] == curr_vertex
                                );
                            }
                        }

                        if new_face.num_vertices_or_edges > 2 {
                            // Sometimes degenerate faces may be generated
                            // due to numerical errors resulting in an isolated
                            // edge not being deleted.
                            //
                            // This kind of degenerate faces are not valid.
                            faces.push(new_face);
                        }
                        break;
                    }
                }
            }
        }

        // Update face ids inside edges so that they point to the faces instead of the triangles.
        for e in &mut edges {
            if let Some(fid) = triangles.get(e.faces[0] as usize)?.parent_face {
                e.faces[0] = fid;
            }

            if let Some(fid) = triangles.get(e.faces[1] as usize)?.parent_face {
                e.faces[1] = fid;
            }
        }

        /*
         * Initialize vertices
         */
        let empty_vertex = Vertex {
            first_adj_face_or_edge: 0,
            num_adj_faces_or_edge: 0,
        };

        vertices.resize(points.len(), empty_vertex);

        // First, find their multiplicities.
        for face in &faces {
            let first_vid = face.first_vertex_or_edge;
            let last_vid = face.first_vertex_or_edge + face.num_vertices_or_edges;

            for i in &vertices_adj_to_face[first_vid as usize..last_vid as usize] {
                vertices[*i as usize].num_adj_faces_or_edge += 1;
            }
        }

        // Now, find their starting id.
        let mut total_num_adj_faces = 0;
        for v in &mut vertices {
            v.first_adj_face_or_edge = total_num_adj_faces;
            total_num_adj_faces += v.num_adj_faces_or_edge;
        }
        faces_adj_to_vertex.resize(total_num_adj_faces as usize, 0);
        edges_adj_to_vertex.resize(total_num_adj_faces as usize, 0);

        // Reset the number of adjascent faces.
        // It will be set againt to the right value as
        // the adjascent face list is filled.
        for v in &mut vertices {
            v.num_adj_faces_or_edge = 0;
        }

        for face_id in 0..faces.len() {
            let face = &faces[face_id];
            let first_vid = face.first_vertex_or_edge;
            let last_vid = face.first_vertex_or_edge + face.num_vertices_or_edges;

            for vid in first_vid..last_vid {
                let v = &mut vertices[vertices_adj_to_face[vid as usize] as usize];
                faces_adj_to_vertex
                    [(v.first_adj_face_or_edge + v.num_adj_faces_or_edge) as usize] =
                    face_id as u32;
                edges_adj_to_vertex
                    [(v.first_adj_face_or_edge + v.num_adj_faces_or_edge) as usize] =
                    edges_adj_to_face[vid as usize];
                v.num_adj_faces_or_edge += 1;
            }
        }

        // Note numerical errors may throw off the Euler characteristic.
        // So we don't check it right now.

        let res = ConvexPolyhedron {
            points,
            vertices,
            faces,
            edges,
            faces_adj_to_vertex,
            edges_adj_to_vertex,
            edges_adj_to_face,
            vertices_adj_to_face,
        };

        // FIXME: for debug.
        // res.check_geometry();

        Some(res)
    }

    /// Verify if this convex polyhedron is actually convex.
    #[inline]
    pub fn check_geometry(&self) {
        for face in &self.faces {
            let p0 =
                self.points[self.vertices_adj_to_face[face.first_vertex_or_edge as usize] as usize];

            for v in &self.points {
                assert!((v - p0).dot(face.normal.as_ref()) <= crate::math::DEFAULT_EPSILON);
            }
        }
    }

    /// The set of vertices of this convex polyhedron.
    #[inline]
    pub fn points(&self) -> &[Point<Real>] {
        &self.points[..]
    }

    /// The topology of the vertices of this convex polyhedron.
    #[inline]
    pub fn vertices(&self) -> &[Vertex] {
        &self.vertices[..]
    }

    /// The topology of the edges of this convex polyhedron.
    #[inline]
    pub fn edges(&self) -> &[Edge] {
        &self.edges[..]
    }

    /// The topology of the faces of this convex polyhedron.
    #[inline]
    pub fn faces(&self) -> &[Face] {
        &self.faces[..]
    }

    /// The array containing the indices of the vertices adjacent to each face.
    #[inline]
    pub fn vertices_adj_to_face(&self) -> &[u32] {
        &self.vertices_adj_to_face[..]
    }

    /// The array containing the indices of the edges adjacent to each face.
    #[inline]
    pub fn edges_adj_to_face(&self) -> &[u32] {
        &self.edges_adj_to_face[..]
    }

    /// The array containing the indices of the faces adjacent to each vertex.
    #[inline]
    pub fn faces_adj_to_vertex(&self) -> &[u32] {
        &self.faces_adj_to_vertex[..]
    }

    fn support_feature_id_toward_eps(
        &self,
        local_dir: &Unit<Vector<Real>>,
        eps: Real,
    ) -> FeatureId {
        let (seps, ceps) = ComplexField::sin_cos(eps);
        let support_pt_id = utils::point_cloud_support_point_id(local_dir.as_ref(), &self.points);
        let vertex = &self.vertices[support_pt_id];

        // Check faces.
        for i in 0..vertex.num_adj_faces_or_edge {
            let face_id = self.faces_adj_to_vertex[(vertex.first_adj_face_or_edge + i) as usize];
            let face = &self.faces[face_id as usize];

            if face.normal.dot(local_dir.as_ref()) >= ceps {
                return FeatureId::Face(face_id);
            }
        }

        // Check edges.
        for i in 0..vertex.num_adj_faces_or_edge {
            let edge_id = self.edges_adj_to_vertex[(vertex.first_adj_face_or_edge + i) as usize];
            let edge = &self.edges[edge_id as usize];

            if edge.dir.dot(local_dir.as_ref()).abs() <= seps {
                return FeatureId::Edge(edge_id);
            }
        }

        // The vertex is the support feature.
        FeatureId::Vertex(support_pt_id as u32)
    }

    /// Computes the ID of the features with a normal that maximize the dot-product with `local_dir`.
    pub fn support_feature_id_toward(&self, local_dir: &Unit<Vector<Real>>) -> FeatureId {
        let eps: Real = na::convert::<f64, Real>(f64::consts::PI / 180.0);
        self.support_feature_id_toward_eps(local_dir, eps)
    }

    /// The normal of the given feature.
    pub fn feature_normal(&self, feature: FeatureId) -> Option<Unit<Vector<Real>>> {
        match feature {
            FeatureId::Face(id) => Some(self.faces[id as usize].normal),
            FeatureId::Edge(id) => {
                let edge = &self.edges[id as usize];
                Some(Unit::new_normalize(
                    *self.faces[edge.faces[0] as usize].normal
                        + *self.faces[edge.faces[1] as usize].normal,
                ))
            }
            FeatureId::Vertex(id) => {
                let vertex = &self.vertices[id as usize];
                let first = vertex.first_adj_face_or_edge;
                let last = vertex.first_adj_face_or_edge + vertex.num_adj_faces_or_edge;
                let mut normal = Vector::zeros();

                for face in &self.faces_adj_to_vertex[first as usize..last as usize] {
                    normal += *self.faces[*face as usize].normal
                }

                Some(Unit::new_normalize(normal))
            }
            FeatureId::Unknown => None,
        }
    }
}

impl SupportMap for ConvexPolyhedron {
    #[inline]
    fn local_support_point(&self, dir: &Vector<Real>) -> Point<Real> {
        utils::point_cloud_support_point(dir, self.points())
    }
}

impl PolygonalFeatureMap for ConvexPolyhedron {
    fn local_support_feature(&self, dir: &Unit<Vector<Real>>, out_feature: &mut PolygonalFeature) {
        let mut best_fid = 0;
        let mut best_dot = self.faces[0].normal.dot(dir);

        for (fid, face) in self.faces[1..].iter().enumerate() {
            let new_dot = face.normal.dot(dir);

            if new_dot > best_dot {
                best_fid = fid + 1;
                best_dot = new_dot;
            }
        }

        let face = &self.faces[best_fid];
        let i1 = face.first_vertex_or_edge;
        // TODO: if there are more than 4 vertices, we need to select four vertices that maximize the area.
        let num_vertices = face.num_vertices_or_edges.min(4);
        let i2 = i1 + num_vertices;

        for (i, (vid, eid)) in self.vertices_adj_to_face[i1 as usize..i2 as usize]
            .iter()
            .zip(self.edges_adj_to_face[i1 as usize..i2 as usize].iter())
            .enumerate()
        {
            out_feature.vertices[i] = self.points[*vid as usize];
            out_feature.vids[i] = *vid;
            out_feature.eids[i] = *eid;
        }

        out_feature.fid = best_fid as u32;
        out_feature.num_vertices = num_vertices as usize;
    }
}

/*
impl ConvexPolyhedron for ConvexPolyhedron {
    fn vertex(&self, id: FeatureId) -> Point<Real> {
        self.points[id.unwrap_vertex() as usize]
    }

    fn edge(&self, id: FeatureId) -> (Point<Real>, Point<Real>, FeatureId, FeatureId) {
        let edge = &self.edges[id.unwrap_edge() as usize];
        let v1 = edge.vertices[0];
        let v2 = edge.vertices[1];

        (
            self.points[v1 as usize],
            self.points[v2 as usize],
            FeatureId::Vertex(v1),
            FeatureId::Vertex(v2),
        )
    }

    fn face(&self, id: FeatureId, out: &mut ConvexPolygonalFeature) {
        out.clear();

        let face = &self.faces[id.unwrap_face() as usize];
        let first_vertex = face.first_vertex_or_edge;
        let last_vertex = face.first_vertex_or_edge + face.num_vertices_or_edges;

        for i in first_vertex..last_vertex {
            let vid = self.vertices_adj_to_face[i];
            let eid = self.edges_adj_to_face[i];
            out.push(self.points[vid], FeatureId::Vertex(vid));
            out.push_edge_feature_id(FeatureId::Edge(eid));
        }

        out.set_normal(face.normal);
        out.set_feature_id(id);
        out.recompute_edge_normals();
    }

    fn support_face_toward(
        &self,
        m: &Isometry<Real>,
        dir: &Unit<Vector<Real>>,
        out: &mut ConvexPolygonalFeature,
    ) {
        let ls_dir = m.inverse_transform_vector(dir);
        let mut best_face = 0;
        let mut max_dot = self.faces[0].normal.dot(&ls_dir);

        for i in 1..self.faces.len() {
            let face = &self.faces[i];
            let dot = face.normal.dot(&ls_dir);

            if dot > max_dot {
                max_dot = dot;
                best_face = i;
            }
        }

        self.face(FeatureId::Face(best_face), out);
        out.transform_by(m);
    }

    fn support_feature_toward(
        &self,
        transform: &Isometry<Real>,
        dir: &Unit<Vector<Real>>,
        angle: Real,
        out: &mut ConvexPolygonalFeature,
    ) {
        out.clear();
        let local_dir = transform.inverse_transform_unit_vector(dir);
        let fid = self.support_feature_id_toward_eps(&local_dir, angle);

        match fid {
            FeatureId::Vertex(_) => {
                let v = self.vertex(fid);
                out.push(v, fid);
                out.set_feature_id(fid);
            }
            FeatureId::Edge(_) => {
                let edge = self.edge(fid);
                out.push(edge.0, edge.2);
                out.push(edge.1, edge.3);
                out.set_feature_id(fid);
                out.push_edge_feature_id(fid);
            }
            FeatureId::Face(_) => self.face(fid, out),
            FeatureId::Unknown => unreachable!(),
        }

        out.transform_by(transform);
    }
}
*/
