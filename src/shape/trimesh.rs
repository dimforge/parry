use std::collections::HashSet;

use crate::bounding_volume::AABB;
use crate::math::{Isometry, Point, Real, Vector};
use crate::partitioning::QBVH;
use crate::shape::composite_shape::SimdCompositeShape;
use crate::shape::{FeatureId, Shape, Triangle, TypedSimdCompositeShape};
#[cfg(feature = "dim2")]
use crate::transformation::ear_clipping::triangulate_ear_clipping;
use crate::utils::hashmap::{Entry, HashMap};
use crate::utils::HashablePartialEq;
#[cfg(feature = "dim3")]
use {crate::shape::Cuboid, crate::utils::SortedPair};

/// Indicated an inconsistency in the topology of a triangle mesh.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum TopologyError {
    /// Found a triangle with two or three identical vertices.
    BadTriangle(u32),
    /// At least two adjacent triangles have opposite orientations.
    BadAdjascentTrianglesOrientation {
        /// The first triangle, with an orientation opposite to the second triangle.
        triangle1: u32,
        /// The second triangle, with an orientation opposite to the first triangle.
        triangle2: u32,
        /// The edge shared between the two triangles.
        edge: (u32, u32),
    },
}

impl std::fmt::Display for TopologyError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::BadTriangle(fid) => {
                f.pad(&format!("the triangle {fid} has at least two identical vertices."))
            }
            Self::BadAdjascentTrianglesOrientation {
                triangle1,
                triangle2,
                edge,
            } => f.pad(&format!("the triangles {triangle1} and {triangle2} sharing the edge {:?} have opposite orientations.", edge)),
        }
    }
}

#[cfg(feature = "std")]
impl std::error::Error for TopologyError {}

/// The set of pseudo-normals of a triangle mesh.
///
/// These pseudo-normals are used for the inside-outside test of a
/// point on the triangle, as described in the paper:
/// "Signed distance computation using the angle weighted pseudonormal", Baerentzen, et al.
/// DOI: 10.1109/TVCG.2005.49
#[derive(Clone)]
#[cfg_attr(feature = "serde-serialize", derive(Serialize, Deserialize))]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Deserialize, rkyv::Serialize)
)]
#[cfg(feature = "dim3")]
pub struct TriMeshPseudoNormals {
    /// The pseudo-normals of the vertices.
    pub vertices_pseudo_normal: Vec<Vector<Real>>,
    /// The pseudo-normals of the edges.
    pub edges_pseudo_normal: HashMap<SortedPair<u32>, Vector<Real>>,
}

#[cfg(feature = "dim3")]
impl Default for TriMeshPseudoNormals {
    fn default() -> Self {
        Self {
            vertices_pseudo_normal: vec![],
            edges_pseudo_normal: Default::default(),
        }
    }
}

/// The connected-components of a triangle mesh.
#[derive(Clone, Debug)]
#[cfg_attr(feature = "serde-serialize", derive(Serialize, Deserialize))]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Deserialize, rkyv::Serialize)
)]
pub struct TriMeshConnectedComponents {
    /// The `face_colors[i]` gives the connected-component index
    /// of the i-th face.
    pub face_colors: Vec<u32>,
    /// The set of faces grouped by connected components.
    pub grouped_faces: Vec<u32>,
    /// The range of connected components. `self.grouped_faces[self.ranges[i]..self.ranges[i + 1]]`
    /// contains the indices of all the faces part of the i-th connected component.
    pub ranges: Vec<usize>,
}

impl TriMeshConnectedComponents {
    /// The total number of connected components.
    pub fn num_connected_components(&self) -> usize {
        self.ranges.len() - 1
    }
}

/// A vertex of a triangle-mesh’s half-edge topology.
#[derive(Clone, Copy, Debug)]
#[cfg_attr(feature = "serde-serialize", derive(Serialize, Deserialize))]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Deserialize, rkyv::Serialize)
)]
pub struct TopoVertex {
    /// One of the half-edge with this vertex as endpoint.
    pub half_edge: u32,
}

/// A face of a triangle-mesh’s half-edge topology.
#[derive(Clone, Copy, Debug)]
#[cfg_attr(feature = "serde-serialize", derive(Serialize, Deserialize))]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Deserialize, rkyv::Serialize)
)]
pub struct TopoFace {
    /// The half-edge adjascent to this face, whith a starting point equal
    /// to the first point of this face.
    pub half_edge: u32,
}

/// A half-edge of a triangle-mesh’s half-edge topology.
#[derive(Clone, Copy, Debug)]
#[cfg_attr(feature = "serde-serialize", derive(Serialize, Deserialize))]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Deserialize, rkyv::Serialize)
)]
pub struct TopoHalfEdge {
    /// The next half-edge.
    pub next: u32,
    /// This half-edge twin on the adjascent triangle.
    ///
    /// This is `u32::MAX` if there is no twin.
    pub twin: u32,
    /// The first vertex of this edge.
    pub vertex: u32,
    /// The face associated to this half-edge.
    pub face: u32,
}

/// The half-edge topology information of a triangle mesh.
#[derive(Clone, Default)]
#[cfg_attr(feature = "serde-serialize", derive(Serialize, Deserialize))]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Deserialize, rkyv::Serialize)
)]
pub struct TriMeshTopology {
    /// The vertices of this half-edge representation.
    pub vertices: Vec<TopoVertex>,
    /// The faces of this half-edge representation.
    pub faces: Vec<TopoFace>,
    /// The half-edges of this half-edge representation.
    pub half_edges: Vec<TopoHalfEdge>,
}

impl TriMeshTopology {
    #[cfg(feature = "dim3")]
    pub(crate) fn face_half_edges_ids(&self, fid: u32) -> [u32; 3] {
        let first_half_edge = self.faces[fid as usize].half_edge;

        let mut result = [first_half_edge; 3];
        for k in 1..3 {
            let half_edge = self.half_edges[result[k - 1] as usize];
            result[k] = half_edge.next;
        }

        result
    }
}

bitflags::bitflags! {
    #[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
    #[cfg_attr(
        feature = "rkyv",
        derive(rkyv::Archive, rkyv::Deserialize, rkyv::Serialize)
    )]
    #[cfg_attr(feature = "cuda", derive(cust_core::DeviceCopy))]
    #[derive(Default)]
    /// The status of the cell of an heightfield.
    pub struct TriMeshFlags: u8 {
        /// If set, the half-edge topology of the trimesh will be computed if possible.
        const HALF_EDGE_TOPOLOGY = 0b0000_0001;
        /// If set, the half-edge topology and connected components of the trimesh will be computed if possible.
        ///
        /// Because of the way it is currently implemented, connected components can only be computed on
        /// a mesh where the half-edge topology computation succeeds. It will no longer be the case in the
        /// future once we decouple the computations.
        const CONNECTED_COMPONENTS = 0b0000_0010;
        /// If set, any triangle that results in a failing half-hedge topology computation will be deleted.
        const DELETE_BAD_TOPOLOGY_TRIANGLES = 0b0000_0100;
        /// If set, the trimesh will be assumed to be oriented (with outward normals).
        ///
        /// The pseudo-normals of its vertices and edges will be computed.
        const ORIENTED = 0b0000_1000;
        /// If set, the duplicate vertices of the trimesh will be merged.
        ///
        /// Two vertices with the exact same coordinates will share the same entry on the
        /// vertex buffer and the index buffer is adjusted accordingly.
        const MERGE_DUPLICATE_VERTICES = 0b0001_0000;
        /// If set, the triangles sharing two vertices with identical index values will be removed.
        ///
        /// Because of the way it is currently implemented, this methods implies that duplicate
        /// vertices will be merged. It will no longer be the case in the future once we decouple
        /// the computations.
        const DELETE_DEGENERATE_TRIANGLES = 0b0010_0000;
        /// If set, two triangles sharing three vertices with identical index values (in any order) will be removed.
        ///
        /// Because of the way it is currently implemented, this methods implies that duplicate
        /// vertices will be merged. It will no longer be the case in the future once we decouple
        /// the computations.
        const DELETE_DUPLICATE_TRIANGLES = 0b0100_0000;
    }
}

#[derive(Clone)]
#[cfg_attr(feature = "serde-serialize", derive(Serialize, Deserialize))]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Deserialize, rkyv::Serialize)
)]
/// A triangle mesh.
pub struct TriMesh {
    qbvh: QBVH<u32>,
    vertices: Vec<Point<Real>>,
    indices: Vec<[u32; 3]>,
    #[cfg(feature = "dim3")]
    pub(crate) pseudo_normals: Option<TriMeshPseudoNormals>,
    topology: Option<TriMeshTopology>,
    connected_components: Option<TriMeshConnectedComponents>,
    flags: TriMeshFlags,
}

impl TriMesh {
    /// Creates a new triangle mesh from a vertex buffer and an index buffer.
    pub fn new(vertices: Vec<Point<Real>>, indices: Vec<[u32; 3]>) -> Self {
        Self::with_flags(vertices, indices, TriMeshFlags::empty())
    }

    /// Creates a new triangle mesh from a vertex buffer and an index buffer, and flags controlling optional properties.
    pub fn with_flags(
        vertices: Vec<Point<Real>>,
        indices: Vec<[u32; 3]>,
        flags: TriMeshFlags,
    ) -> Self {
        assert!(
            indices.len() > 0,
            "A triangle mesh must contain at least one triangle."
        );

        let mut result = Self {
            qbvh: QBVH::new(),
            vertices,
            indices,
            #[cfg(feature = "dim3")]
            pseudo_normals: None,
            topology: None,
            connected_components: None,
            flags: TriMeshFlags::empty(),
        };

        let _ = result.set_flags(flags);

        if result.qbvh.raw_nodes().is_empty() {
            // The QBVH hasn’t been computed by `.set_flags`.
            result.rebuild_qbvh();
        }

        result
    }

    /// The flags of this triangle mesh.
    pub fn flags(&self) -> TriMeshFlags {
        self.flags
    }

    /// Sets the flags of this triangle mesh, controlling its optional associated data.
    pub fn set_flags(&mut self, flags: TriMeshFlags) -> Result<(), TopologyError> {
        let mut result = Ok(());
        let prev_indices_len = self.indices.len();

        if !flags.contains(TriMeshFlags::HALF_EDGE_TOPOLOGY) {
            self.topology = None;
        }

        #[cfg(feature = "dim3")]
        if !flags.contains(TriMeshFlags::ORIENTED) {
            self.pseudo_normals = None;
        }

        if !flags.contains(TriMeshFlags::CONNECTED_COMPONENTS) {
            self.connected_components = None;
        }

        let difference = flags & !self.flags;

        if difference.intersects(
            TriMeshFlags::MERGE_DUPLICATE_VERTICES
                | TriMeshFlags::DELETE_DEGENERATE_TRIANGLES
                | TriMeshFlags::DELETE_DUPLICATE_TRIANGLES,
        ) {
            self.merge_duplicate_vertices(
                flags.contains(TriMeshFlags::DELETE_DEGENERATE_TRIANGLES),
                flags.contains(TriMeshFlags::DELETE_DUPLICATE_TRIANGLES),
            )
        }

        if difference.intersects(
            TriMeshFlags::HALF_EDGE_TOPOLOGY
                | TriMeshFlags::CONNECTED_COMPONENTS
                | TriMeshFlags::DELETE_BAD_TOPOLOGY_TRIANGLES,
        ) {
            result = self.compute_topology(
                flags.contains(TriMeshFlags::CONNECTED_COMPONENTS),
                flags.contains(TriMeshFlags::DELETE_BAD_TOPOLOGY_TRIANGLES),
            );
        }

        #[cfg(feature = "dim3")]
        if difference.contains(TriMeshFlags::ORIENTED) {
            self.compute_pseudo_normals();
        }

        if prev_indices_len != self.indices.len() {
            self.rebuild_qbvh();
        }

        self.flags = flags;
        result
    }

    /// Transforms in-place the vertices of this triangle mesh.
    pub fn transform_vertices(&mut self, transform: &Isometry<Real>) {
        self.vertices
            .iter_mut()
            .for_each(|pt| *pt = transform * *pt);
        self.rebuild_qbvh();

        // The pseudo-normals must be rotated too.
        #[cfg(feature = "dim3")]
        if let Some(pseudo_normals) = &mut self.pseudo_normals {
            pseudo_normals
                .vertices_pseudo_normal
                .iter_mut()
                .for_each(|n| *n = transform * *n);
            pseudo_normals
                .edges_pseudo_normal
                .values_mut()
                .for_each(|n| *n = transform * *n);
        }
    }

    /// Returns a scaled version of this triangle mesh.
    pub fn scaled(mut self, scale: &Vector<Real>) -> Self {
        self.vertices
            .iter_mut()
            .for_each(|pt| pt.coords.component_mul_assign(scale));

        #[cfg(feature = "dim3")]
        if let Some(pn) = &mut self.pseudo_normals {
            pn.vertices_pseudo_normal.iter_mut().for_each(|n| {
                n.component_mul_assign(scale);
                let _ = n.try_normalize_mut(0.0);
            });
            pn.edges_pseudo_normal.values_mut().for_each(|n| {
                n.component_mul_assign(scale);
                let _ = n.try_normalize_mut(0.0);
            });
        }

        Self {
            qbvh: self.qbvh.scaled(scale),
            vertices: self.vertices,
            indices: self.indices,
            #[cfg(feature = "dim3")]
            pseudo_normals: self.pseudo_normals,
            topology: self.topology,
            connected_components: self.connected_components,
            flags: self.flags,
        }
    }

    /// Appends a second triangle mesh to this triangle mesh.
    pub fn append(&mut self, rhs: &TriMesh) {
        let base_id = self.vertices.len() as u32;
        self.vertices.extend_from_slice(rhs.vertices());
        self.indices.extend(
            rhs.indices()
                .iter()
                .map(|idx| [idx[0] + base_id, idx[1] + base_id, idx[2] + base_id]),
        );

        let vertices = std::mem::replace(&mut self.vertices, Vec::new());
        let indices = std::mem::replace(&mut self.indices, Vec::new());
        *self = TriMesh::with_flags(vertices, indices, self.flags);
    }

    /// Create a `TriMesh` from a set of points assumed to describe a counter-clockwise non-convex polygon.
    ///
    /// This operation may fail if the input polygon is invalid, e.g. it is non-simple or has zero surface area.
    #[cfg(feature = "dim2")]
    pub fn from_polygon(vertices: Vec<Point<Real>>) -> Option<Self> {
        triangulate_ear_clipping(&vertices).map(|indices| Self::new(vertices, indices))
    }

    /// Compute the axis-aligned bounding box of this triangle mesh.
    pub fn aabb(&self, pos: &Isometry<Real>) -> AABB {
        self.qbvh.root_aabb().transform_by(pos)
    }

    /// Gets the local axis-aligned bounding box of this triangle mesh.
    pub fn local_aabb(&self) -> &AABB {
        self.qbvh.root_aabb()
    }

    /// The acceleration structure used by this triangle-mesh.
    pub fn qbvh(&self) -> &QBVH<u32> {
        &self.qbvh
    }

    /// The number of triangles forming this mesh.
    pub fn num_triangles(&self) -> usize {
        self.indices.len()
    }

    /// Does the given feature ID identify a backface of this trimesh?
    pub fn is_backface(&self, feature: FeatureId) -> bool {
        if let FeatureId::Face(i) = feature {
            i >= self.indices.len() as u32
        } else {
            false
        }
    }

    /// An iterator through all the triangles of this mesh.
    pub fn triangles(&self) -> impl ExactSizeIterator<Item = Triangle> + '_ {
        self.indices.iter().map(move |ids| {
            Triangle::new(
                self.vertices[ids[0] as usize],
                self.vertices[ids[1] as usize],
                self.vertices[ids[2] as usize],
            )
        })
    }

    /// Get the `i`-th triangle of this mesh.
    pub fn triangle(&self, i: u32) -> Triangle {
        let idx = self.indices[i as usize];
        Triangle::new(
            self.vertices[idx[0] as usize],
            self.vertices[idx[1] as usize],
            self.vertices[idx[2] as usize],
        )
    }

    /// The vertex buffer of this mesh.
    pub fn vertices(&self) -> &[Point<Real>] {
        &self.vertices[..]
    }

    /// The index buffer of this mesh.
    pub fn indices(&self) -> &[[u32; 3]] {
        &self.indices
    }

    /// A flat view of the index buffer of this mesh.
    pub fn flat_indices(&self) -> &[u32] {
        unsafe {
            let len = self.indices.len() * 3;
            let data = self.indices.as_ptr() as *const u32;
            std::slice::from_raw_parts(data, len)
        }
    }

    /// Returns the topology information of this trimesh, if it has been computed.
    pub fn topology(&self) -> Option<&TriMeshTopology> {
        self.topology.as_ref()
    }

    /// Returns the connected-component information of this trimesh, if it has been computed.
    pub fn connected_components(&self) -> Option<&TriMeshConnectedComponents> {
        self.connected_components.as_ref()
    }

    /// The pseudo-normals of this triangle mesh, if they have been computed.
    #[cfg(feature = "dim3")]
    pub fn pseudo_normals(&self) -> Option<&TriMeshPseudoNormals> {
        self.pseudo_normals.as_ref()
    }

    fn rebuild_qbvh(&mut self) {
        let data = self.indices.iter().enumerate().map(|(i, idx)| {
            let aabb = Triangle::new(
                self.vertices[idx[0] as usize],
                self.vertices[idx[1] as usize],
                self.vertices[idx[2] as usize],
            )
            .local_aabb();
            (i as u32, aabb)
        });

        // NOTE: we apply no dilation factor because we won't
        // update this tree dynamically.
        self.qbvh.clear_and_rebuild(data, 0.0);
    }

    /// Reverse the orientation of the triangle mesh.
    pub fn reverse(&mut self) {
        self.indices.iter_mut().for_each(|idx| idx.swap(0, 1));

        // NOTE: the QBVH, and connected components are not changed by this operation.
        //       The pseudo-normals just have to be flipped.
        //       The topology must be recomputed.

        #[cfg(feature = "dim3")]
        if let Some(pseudo_normals) = &mut self.pseudo_normals {
            for n in &mut pseudo_normals.vertices_pseudo_normal {
                *n = -*n;
            }

            for n in pseudo_normals.edges_pseudo_normal.values_mut() {
                *n = -*n;
            }
        }

        if self.flags.contains(TriMeshFlags::HALF_EDGE_TOPOLOGY) {
            // TODO: this could be done more efficiently.
            let _ = self.compute_topology(false, false);
        }
    }

    /// Merge all duplicate vertices and adjust the index buffer accordingly.
    ///
    /// If `delete_degenerate_triangles` is set to true, any triangle with two
    /// identical vertices will be removed.
    ///
    /// This is typically used to recover a vertex buffer from which we can deduce
    /// adjacency information. between triangles by observing how the vertices are
    /// shared by triangles based on the index buffer.
    fn merge_duplicate_vertices(
        &mut self,
        delete_degenerate_triangles: bool,
        delete_duplicate_triangles: bool,
    ) {
        let mut vtx_to_id = HashMap::default();
        let mut new_vertices = Vec::with_capacity(self.vertices.len());
        let mut new_indices = Vec::with_capacity(self.indices.len());
        let mut triangle_set = HashSet::new();

        fn resolve_coord_id(
            coord: &Point<Real>,
            vtx_to_id: &mut HashMap<HashablePartialEq<Point<Real>>, u32>,
            new_vertices: &mut Vec<Point<Real>>,
        ) -> u32 {
            let key = HashablePartialEq::new(coord.clone());
            let id = match vtx_to_id.entry(key) {
                Entry::Occupied(entry) => entry.into_mut(),
                Entry::Vacant(entry) => entry.insert(new_vertices.len() as u32),
            };

            if *id == new_vertices.len() as u32 {
                new_vertices.push(coord.clone());
            }

            *id
        }

        for t in self.indices.iter() {
            let va = resolve_coord_id(
                &self.vertices[t[0] as usize],
                &mut vtx_to_id,
                &mut new_vertices,
            );

            let vb = resolve_coord_id(
                &self.vertices[t[1] as usize],
                &mut vtx_to_id,
                &mut new_vertices,
            );

            let vc = resolve_coord_id(
                &self.vertices[t[2] as usize],
                &mut vtx_to_id,
                &mut new_vertices,
            );

            let is_degenerate = va == vb || va == vc || vb == vc;

            if !is_degenerate || !delete_degenerate_triangles {
                if delete_duplicate_triangles {
                    let (c, b, a) = crate::utils::sort3(&va, &vb, &vc);
                    if triangle_set.insert((*a, *b, *c)) {
                        new_indices.push([va, vb, vc])
                    }
                } else {
                    new_indices.push([va, vb, vc]);
                }
            }
        }

        new_vertices.shrink_to_fit();

        self.vertices = new_vertices;
        self.indices = new_indices;

        // Vertices and indices changed: the pseudo-normals are no longer valid.
        #[cfg(feature = "dim3")]
        if self.pseudo_normals.is_some() {
            self.compute_pseudo_normals();
        }

        // Vertices and indices changed: the topology no longer valid.
        #[cfg(feature = "dim3")]
        if self.topology.is_some() {
            let _ = self.compute_topology(self.connected_components.is_some(), false);
        }
    }

    #[cfg(feature = "dim3")]
    /// Computes the pseudo-normals used for solid point-projection.
    ///
    /// This computes the pseudo-normals needed by the point containment test described in
    /// "Signed distance computation using the angle weighted pseudonormal", Baerentzen, et al.
    /// DOI: 10.1109/TVCG.2005.49
    ///
    /// For the point-containment test to properly detect the inside of the trimesh (i.e. to return
    /// `proj.is_inside = true`), the trimesh must:
    /// - Be manifold (closed, no t-junctions, etc.)
    /// - Be oriented with outward normals.
    ///
    /// If the the trimesh is correctly oriented, but is manifold everywhere except at its boundaries,
    /// then the computed pseudo-normals will provide correct point-containment test results except
    /// for points closest to the boundary of the mesh.
    ///
    /// It may be useful to call `self.remove_duplicate_vertices()` before this method, in order to fix the
    /// index buffer if some of the vertices of this trimesh are duplicated.
    fn compute_pseudo_normals(&mut self) {
        use na::RealField;

        let mut degenerate_triangles = vec![false; self.indices().len()];
        let mut vertices_pseudo_normal = vec![Vector::zeros(); self.vertices().len()];
        let mut edges_pseudo_normal = HashMap::default();
        let mut edges_multiplicity = HashMap::default();

        for (k, idx) in self.indices().iter().enumerate() {
            let vtx = self.vertices();
            let tri = Triangle::new(
                vtx[idx[0] as usize],
                vtx[idx[1] as usize],
                vtx[idx[2] as usize],
            );

            if let Some(n) = tri.normal() {
                let ang1 = (tri.b - tri.a).angle(&(tri.c - tri.a));
                let ang2 = (tri.a - tri.b).angle(&(tri.c - tri.b));
                let ang3 = (tri.b - tri.c).angle(&(tri.a - tri.c));

                // NOTE: almost-degenerate triangles can result in a bad normal (due to
                //       float errors, with a high weight (due to the large angle of one
                //       of its vertices). We need to ignore these triangles.
                degenerate_triangles[k] = ang1.max(ang2).max(ang3) > Real::pi() - Real::pi() / 10.0;

                if degenerate_triangles[k] {
                    continue;
                }

                vertices_pseudo_normal[idx[0] as usize] += *n * ang1;
                vertices_pseudo_normal[idx[1] as usize] += *n * ang2;
                vertices_pseudo_normal[idx[2] as usize] += *n * ang3;

                let edges = [
                    SortedPair::new(idx[0], idx[1]),
                    SortedPair::new(idx[0], idx[2]),
                    SortedPair::new(idx[1], idx[2]),
                ];

                for edge in &edges {
                    let edge_n = edges_pseudo_normal
                        .entry(*edge)
                        .or_insert_with(Vector::zeros);
                    *edge_n += *n; // NOTE: there is no need to multiply by the incident angle since it is always equal to PI for all the edges.
                    let edge_mult = edges_multiplicity.entry(*edge).or_insert(0);
                    *edge_mult += 1;
                }
            }
        }

        self.pseudo_normals = Some(TriMeshPseudoNormals {
            vertices_pseudo_normal,
            edges_pseudo_normal,
        })
    }

    fn delete_bad_topology_triangles(&mut self) {
        let mut half_edge_set = HashSet::new();
        let mut deleted_any = false;

        // First, create three half-edges for each face.
        self.indices.retain(|idx| {
            if idx[0] == idx[1] || idx[0] == idx[2] || idx[1] == idx[2] {
                deleted_any = true;
                return false;
            }

            for k in 0..3 {
                let edge_key = (idx[k as usize], idx[(k as usize + 1) % 3]);
                if half_edge_set.contains(&edge_key) {
                    deleted_any = true;
                    return false;
                }
            }

            for k in 0..3 {
                let edge_key = (idx[k as usize], idx[(k as usize + 1) % 3]);
                let _ = half_edge_set.insert(edge_key);
            }

            true
        });
    }

    /// Computes half-edge topological information for this triangle mesh, based on its index buffer only.
    ///
    /// This computes the half-edge representation of this triangle mesh’s topology. This is useful for advanced
    /// geometric operations like trimesh-trimesh intersection geometry computation.
    ///
    /// It may be useful to call `self.merge_duplicate_vertices(true, true)` before this method, in order to fix the
    /// index buffer if some of the vertices of this trimesh are duplicated.
    ///
    /// # Return
    /// Returns `true` if the computation succeeded. Returns `false` if this mesh can’t have an half-edge representation
    /// because at least three faces share the same edge.
    fn compute_topology(
        &mut self,
        compute_connected_components: bool,
        delete_bad_triangles: bool,
    ) -> Result<(), TopologyError> {
        if delete_bad_triangles {
            self.delete_bad_topology_triangles();
        }

        let mut topology = TriMeshTopology::default();
        let mut half_edge_map = HashMap::default();

        topology.vertices.resize(
            self.vertices.len(),
            TopoVertex {
                half_edge: u32::MAX,
            },
        );

        // First, create three half-edges for each face.
        for (fid, idx) in self.indices.iter().enumerate() {
            let half_edge_base_id = topology.half_edges.len() as u32;

            if idx[0] == idx[1] || idx[0] == idx[2] || idx[1] == idx[2] {
                return Err(TopologyError::BadTriangle(fid as u32));
            }

            for k in 0u32..3 {
                let half_edge = TopoHalfEdge {
                    next: half_edge_base_id + (k + 1) % 3,
                    // We don’t know which one it is yet.
                    // If the twin doesn’t exist, we use `u32::MAX` as
                    // it’s (invalid) index. This value can be relied on
                    // by other algorithms.
                    twin: u32::MAX,
                    vertex: idx[k as usize],
                    face: fid as u32,
                };
                topology.half_edges.push(half_edge);

                let edge_key = (idx[k as usize], idx[(k as usize + 1) % 3]);

                if let Some(existing) = half_edge_map.insert(edge_key, half_edge_base_id + k) {
                    // If the same edge already exists (with the same vertex order) then
                    // we have two triangles sharing the same but with opposite incompatible orientations.
                    return Err(TopologyError::BadAdjascentTrianglesOrientation {
                        edge: edge_key,
                        triangle1: topology.half_edges[existing as usize].face,
                        triangle2: fid as u32,
                    });
                }

                topology.vertices[idx[k as usize] as usize].half_edge = half_edge_base_id + k;
            }

            topology.faces.push(TopoFace {
                half_edge: half_edge_base_id,
            })
        }

        // Second, identify twins.
        for (key, he1) in &half_edge_map {
            if key.0 < key.1 {
                // Test, to avoid checking the same pair twice.
                if let Some(he2) = half_edge_map.get(&(key.1, key.0)) {
                    topology.half_edges[*he1 as usize].twin = *he2;
                    topology.half_edges[*he2 as usize].twin = *he1;
                }
            }
        }

        self.topology = Some(topology);

        if compute_connected_components {
            self.compute_connected_components();
        }

        Ok(())
    }

    // NOTE: we can only compute connected components if the topology
    //       has been computed too. So instead of making this method
    //       public, the `.compute_topology` method has a boolean to
    //       compute the connected componets too.
    fn compute_connected_components(&mut self) {
        let topo = self.topology.as_ref().unwrap();
        let mut face_colors = vec![u32::MAX; topo.faces.len()];
        let mut grouped_faces = Vec::new();
        let mut ranges = vec![0];

        for i in 0..topo.faces.len() {
            if face_colors[i] == u32::MAX {
                let mut stack = vec![i as u32];
                let color = ranges.len() as u32 - 1;

                while let Some(tri_id) = stack.pop() {
                    grouped_faces.push(tri_id);
                    let eid = topo.faces[tri_id as usize].half_edge;
                    let edge_a = &topo.half_edges[eid as usize];
                    let edge_b = &topo.half_edges[edge_a.next as usize];
                    let edge_c = &topo.half_edges[edge_b.next as usize];
                    let edges = [edge_a, edge_b, edge_c];

                    for edge in edges {
                        if let Some(twin) = topo.half_edges.get(edge.twin as usize) {
                            if face_colors[twin.face as usize] == u32::MAX {
                                face_colors[twin.face as usize] = color;
                                grouped_faces.push(twin.face);
                                stack.push(twin.face);
                            }
                        }
                    }
                }

                ranges.push(grouped_faces.len());
            }
        }

        self.connected_components = Some(TriMeshConnectedComponents {
            face_colors,
            grouped_faces,
            ranges,
        });
    }

    #[allow(dead_code)] // Useful for testing.
    pub(crate) fn assert_half_edge_topology_is_valid(&self) {
        let topo = self
            .topology
            .as_ref()
            .expect("No topology information found.");
        assert_eq!(self.vertices.len(), topo.vertices.len());
        assert_eq!(self.indices.len(), topo.faces.len());

        for (face_id, (face, idx)) in topo.faces.iter().zip(self.indices.iter()).enumerate() {
            let he0 = topo.half_edges[face.half_edge as usize];
            assert_eq!(he0.face, face_id as u32);
            assert_eq!(he0.vertex, idx[0]);
            let he1 = topo.half_edges[he0.next as usize];
            assert_eq!(he1.face, face_id as u32);
            assert_eq!(he1.vertex, idx[1]);
            let he2 = topo.half_edges[he1.next as usize];
            assert_eq!(he2.face, face_id as u32);
            assert_eq!(he2.vertex, idx[2]);
            assert_eq!(he2.next, face.half_edge);
        }

        for he in &topo.half_edges {
            let idx = &self.indices[he.face as usize];
            assert!(he.vertex == idx[0] || he.vertex == idx[1] || he.vertex == idx[2]);
        }
    }
}

/*
#[cfg(feature = "dim3")]
impl RayCast for TriMesh {
    fn cast_local_ray_and_get_normal(
        &self,
        ray: &Ray,
        max_toi: Real,
        solid: bool,
    ) -> Option<RayIntersection> {
        // FIXME: do a best-first search.
        let mut intersections = Vec::new();
        self.qbvh.cast_ray(&ray, max_toi, &mut intersections);
        let mut best: Option<RayIntersection> = None;

        for inter in intersections {
            let tri = self.triangle(inter);
            if let Some(inter) = tri.cast_local_ray_and_get_normal(ray, max_toi, solid) {
                if let Some(curr) = &mut best {
                    if curr.toi > inter.toi {
                        *curr = inter;
                    }
                } else {
                    best = Some(inter);
                }
            }
        }

        best
    }

    fn intersects_local_ray(&self, ray: &Ray, max_toi: Real) -> bool {
        // FIXME: do a best-first search.
        let mut intersections = Vec::new();
        self.qbvh.cast_ray(&ray, max_toi, &mut intersections);

        for inter in intersections {
            let tri = self.triangle(inter);
            if tri.intersects_local_ray(ray, max_toi) {
                return true;
            }
        }

        false
    }
}
*/

#[cfg(feature = "dim3")]
impl<Heights, Status> From<crate::shape::GenericHeightField<Heights, Status>> for TriMesh
where
    Heights: crate::shape::HeightFieldStorage<Item = Real>,
    Status: crate::shape::HeightFieldStorage<Item = crate::shape::HeightFieldCellStatus>,
{
    fn from(heightfield: crate::shape::GenericHeightField<Heights, Status>) -> Self {
        let (vtx, idx) = heightfield.to_trimesh();
        TriMesh::new(vtx, idx)
    }
}

#[cfg(feature = "dim3")]
impl From<Cuboid> for TriMesh {
    fn from(cuboid: Cuboid) -> Self {
        let (vtx, idx) = cuboid.to_trimesh();
        TriMesh::new(vtx, idx)
    }
}

impl SimdCompositeShape for TriMesh {
    fn map_part_at(&self, i: u32, f: &mut dyn FnMut(Option<&Isometry<Real>>, &dyn Shape)) {
        let tri = self.triangle(i);
        f(None, &tri)
    }

    fn qbvh(&self) -> &QBVH<u32> {
        &self.qbvh
    }
}

impl TypedSimdCompositeShape for TriMesh {
    type PartShape = Triangle;
    type PartId = u32;

    #[inline(always)]
    fn map_typed_part_at(
        &self,
        i: u32,
        mut f: impl FnMut(Option<&Isometry<Real>>, &Self::PartShape),
    ) {
        let tri = self.triangle(i);
        f(None, &tri)
    }

    #[inline(always)]
    fn map_untyped_part_at(&self, i: u32, mut f: impl FnMut(Option<&Isometry<Real>>, &dyn Shape)) {
        let tri = self.triangle(i);
        f(None, &tri)
    }

    fn typed_qbvh(&self) -> &QBVH<u32> {
        &self.qbvh
    }
}
