use crate::bounding_volume::Aabb;
use crate::math::{Isometry, Point, Real, Vector};
use crate::partitioning::Qbvh;
use crate::shape::{FeatureId, Shape, Triangle, TrianglePseudoNormals, TypedSimdCompositeShape};
use crate::utils::HashablePartialEq;
use alloc::{vec, vec::Vec};
use core::fmt;
#[cfg(feature = "dim3")]
use {crate::shape::Cuboid, crate::utils::SortedPair, na::Unit};

use {
    crate::shape::composite_shape::SimdCompositeShape,
    crate::utils::hashmap::{Entry, HashMap},
    crate::utils::hashset::HashSet,
};

#[cfg(feature = "dim2")]
use crate::transformation::ear_clipping::triangulate_ear_clipping;

use crate::query::details::NormalConstraints;
#[cfg(feature = "rkyv")]
use rkyv::{bytecheck, CheckBytes};

/// Indicated an inconsistency in the topology of a triangle mesh.
#[derive(thiserror::Error, Copy, Clone, Debug, PartialEq, Eq)]
pub enum TopologyError {
    /// Found a triangle with two or three identical vertices.
    #[error("the triangle {0} has at least two identical vertices.")]
    BadTriangle(u32),
    /// At least two adjacent triangles have opposite orientations.
    #[error("the triangles {triangle1} and {triangle2} sharing the edge {edge:?} have opposite orientations.")]
    BadAdjacentTrianglesOrientation {
        /// The first triangle, with an orientation opposite to the second triangle.
        triangle1: u32,
        /// The second triangle, with an orientation opposite to the first triangle.
        triangle2: u32,
        /// The edge shared between the two triangles.
        edge: (u32, u32),
    },
}

/// Indicated an inconsistency while building a triangle mesh.
#[derive(thiserror::Error, Copy, Clone, Debug, PartialEq, Eq)]
pub enum TriMeshBuilderError {
    /// A triangle mesh must contain at least one triangle.
    #[error("A triangle mesh must contain at least one triangle.")]
    EmptyIndices,
    /// Indicated an inconsistency in the topology of a triangle mesh.
    #[error("Topology Error: {0}")]
    TopologyError(TopologyError),
}

/// The set of pseudo-normals of a triangle mesh.
///
/// These pseudo-normals are used for the inside-outside test of a
/// point on the triangle, as described in the paper:
/// "Signed distance computation using the angle weighted pseudonormal", Baerentzen, et al.
/// DOI: 10.1109/TVCG.2005.49
#[derive(Default, Clone)]
#[cfg_attr(feature = "serde-serialize", derive(Serialize, Deserialize))]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Deserialize, rkyv::Serialize),
    archive(check_bytes)
)]
#[repr(C)]
#[cfg(feature = "dim3")]
pub struct TriMeshPseudoNormals {
    /// The pseudo-normals of the vertices.
    pub vertices_pseudo_normal: Vec<Vector<Real>>,
    /// The pseudo-normals of the edges.
    pub edges_pseudo_normal: Vec<[Vector<Real>; 3]>,
}

/// The connected-components of a triangle mesh.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde-serialize", derive(Serialize, Deserialize))]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Deserialize, rkyv::Serialize),
    archive(check_bytes)
)]
#[repr(C)]
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

    /// Convert the connected-component description into actual meshes (returned as raw index and
    /// vertex buffers).
    ///
    /// The `mesh` must be the one used to generate `self`, otherwise it might panic or produce an
    /// unexpected result.
    pub fn to_mesh_buffers(&self, mesh: &TriMesh) -> Vec<(Vec<Point<Real>>, Vec<[u32; 3]>)> {
        let mut result = vec![];
        let mut new_vtx_index: Vec<_> = vec![u32::MAX; mesh.vertices.len()];

        for ranges in self.ranges.windows(2) {
            let num_faces = ranges[1] - ranges[0];

            if num_faces == 0 {
                continue;
            }

            let mut vertices = Vec::with_capacity(num_faces);
            let mut indices = Vec::with_capacity(num_faces);

            for fid in ranges[0]..ranges[1] {
                let vids = mesh.indices[self.grouped_faces[fid] as usize];
                let new_vids = vids.map(|id| {
                    if new_vtx_index[id as usize] == u32::MAX {
                        vertices.push(mesh.vertices[id as usize]);
                        new_vtx_index[id as usize] = vertices.len() as u32 - 1;
                    }

                    new_vtx_index[id as usize]
                });
                indices.push(new_vids);
            }

            result.push((vertices, indices));
        }

        result
    }

    /// Convert the connected-component description into actual meshes.
    ///
    /// The `mesh` must be the one used to generate `self`, otherwise it might panic or produce an
    /// unexpected result.
    ///
    /// All the meshes are constructed with the given `flags`.
    pub fn to_meshes(
        &self,
        mesh: &TriMesh,
        flags: TriMeshFlags,
    ) -> Vec<Result<TriMesh, TriMeshBuilderError>> {
        self.to_mesh_buffers(mesh)
            .into_iter()
            .map(|(vtx, idx)| TriMesh::with_flags(vtx, idx, flags))
            .collect()
    }
}

/// A vertex of a triangle-mesh’s half-edge topology.
#[derive(Clone, Copy, Debug)]
#[cfg_attr(feature = "serde-serialize", derive(Serialize, Deserialize))]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Deserialize, rkyv::Serialize, CheckBytes),
    archive(as = "Self")
)]
#[repr(C)]
pub struct TopoVertex {
    /// One of the half-edge with this vertex as endpoint.
    pub half_edge: u32,
}

/// A face of a triangle-mesh’s half-edge topology.
#[derive(Clone, Copy, Debug)]
#[cfg_attr(feature = "serde-serialize", derive(Serialize, Deserialize))]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Deserialize, rkyv::Serialize, CheckBytes),
    archive(as = "Self")
)]
#[repr(C)]
pub struct TopoFace {
    /// The half-edge adjacent to this face, with a starting point equal
    /// to the first point of this face.
    pub half_edge: u32,
}

/// A half-edge of a triangle-mesh’s half-edge topology.
#[derive(Clone, Copy, Debug)]
#[cfg_attr(feature = "serde-serialize", derive(Serialize, Deserialize))]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Deserialize, rkyv::Serialize, CheckBytes),
    archive(as = "Self")
)]
#[repr(C)]
pub struct TopoHalfEdge {
    /// The next half-edge.
    pub next: u32,
    /// This half-edge twin on the adjacent triangle.
    ///
    /// This is `u32::MAX` if there is no twin.
    pub twin: u32,
    /// The first vertex of this edge.
    pub vertex: u32,
    /// The face associated to this half-edge.
    pub face: u32,
}

/// The half-edge topology information of a triangle mesh.
#[derive(Default, Clone)]
#[cfg_attr(feature = "serde-serialize", derive(Serialize, Deserialize))]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Deserialize, rkyv::Serialize),
    archive(check_bytes)
)]
#[repr(C)]
pub struct TriMeshTopology {
    /// The vertices of this half-edge representation.
    pub vertices: Vec<TopoVertex>,
    /// The faces of this half-edge representation.
    pub faces: Vec<TopoFace>,
    /// The half-edges of this half-edge representation.
    pub half_edges: Vec<TopoHalfEdge>,
}

#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Deserialize, rkyv::Serialize),
    archive(as = "Self")
)]
#[repr(C)]
#[derive(Clone, Copy, Debug, Default, Eq, Hash, Ord, PartialEq, PartialOrd)]
/// Controls how a [`TriMesh`] should be loaded.
pub struct TriMeshFlags(u16);

bitflags::bitflags! {
    impl TriMeshFlags: u16 {
        /// If set, the half-edge topology of the trimesh will be computed if possible.
        const HALF_EDGE_TOPOLOGY = 1;
        /// If set, the connected components of the trimesh will be computed.
        const CONNECTED_COMPONENTS = 1 << 1;
        /// If set, any triangle that results in a failing half-hedge topology computation will be deleted.
        const DELETE_BAD_TOPOLOGY_TRIANGLES = 1 << 2;
        /// If set, the trimesh will be assumed to be oriented (with outward normals).
        ///
        /// The pseudo-normals of its vertices and edges will be computed.
        const ORIENTED = 1 << 3;
        /// If set, the duplicate vertices of the trimesh will be merged.
        ///
        /// Two vertices with the exact same coordinates will share the same entry on the
        /// vertex buffer and the index buffer is adjusted accordingly.
        const MERGE_DUPLICATE_VERTICES = 1 << 4;
        /// If set, the triangles sharing two vertices with identical index values will be removed.
        ///
        /// Because of the way it is currently implemented, this methods implies that duplicate
        /// vertices will be merged. It will no longer be the case in the future once we decouple
        /// the computations.
        const DELETE_DEGENERATE_TRIANGLES = 1 << 5;
        /// If set, two triangles sharing three vertices with identical index values (in any order)
        /// will be removed.
        ///
        /// Because of the way it is currently implemented, this methods implies that duplicate
        /// vertices will be merged. It will no longer be the case in the future once we decouple
        /// the computations.
        const DELETE_DUPLICATE_TRIANGLES = 1 << 6;
        /// If set, a special treatment will be applied to contact manifold calculation to eliminate
        /// or fix contacts normals that could lead to incorrect bumps in physics simulation
        /// (especially on flat surfaces).
        ///
        /// This is achieved by taking into account adjacent triangle normals when computing contact
        /// points for a given triangle.
        const FIX_INTERNAL_EDGES = (1 << 7) | Self::MERGE_DUPLICATE_VERTICES.bits();
    }
}

#[cfg_attr(feature = "serde-serialize", derive(Serialize, Deserialize))]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Deserialize, rkyv::Serialize),
    archive(check_bytes)
)]
#[repr(C)]
#[derive(Clone)]
/// A triangle mesh.
pub struct TriMesh {
    qbvh: Qbvh<u32>,
    vertices: Vec<Point<Real>>,
    indices: Vec<[u32; 3]>,
    #[cfg(feature = "dim3")]
    pub(crate) pseudo_normals: Option<TriMeshPseudoNormals>,
    topology: Option<TriMeshTopology>,
    connected_components: Option<TriMeshConnectedComponents>,
    flags: TriMeshFlags,
}

impl fmt::Debug for TriMesh {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "GenericTriMesh")
    }
}

impl TriMesh {
    /// Creates a new triangle mesh from a vertex buffer and an index buffer.
    pub fn new(
        vertices: Vec<Point<Real>>,
        indices: Vec<[u32; 3]>,
    ) -> Result<Self, TriMeshBuilderError> {
        Self::with_flags(vertices, indices, TriMeshFlags::empty())
    }

    /// Creates a new triangle mesh from a vertex buffer and an index buffer, and flags controlling optional properties.
    pub fn with_flags(
        vertices: Vec<Point<Real>>,
        indices: Vec<[u32; 3]>,
        flags: TriMeshFlags,
    ) -> Result<Self, TriMeshBuilderError> {
        if indices.is_empty() {
            return Err(TriMeshBuilderError::EmptyIndices);
        }

        let mut result = Self {
            qbvh: Qbvh::new(),
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
            // The Qbvh hasn’t been computed by `.set_flags`.
            result.rebuild_qbvh();
        }

        Ok(result)
    }

    /// Sets the flags of this triangle mesh, controlling its optional associated data.
    pub fn set_flags(&mut self, flags: TriMeshFlags) -> Result<(), TopologyError> {
        let mut result = Ok(());
        let prev_indices_len = self.indices.len();

        if !flags.contains(TriMeshFlags::HALF_EDGE_TOPOLOGY) {
            self.topology = None;
        }

        #[cfg(feature = "dim3")]
        if !flags.intersects(TriMeshFlags::ORIENTED | TriMeshFlags::FIX_INTERNAL_EDGES) {
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
            TriMeshFlags::HALF_EDGE_TOPOLOGY | TriMeshFlags::DELETE_BAD_TOPOLOGY_TRIANGLES,
        ) {
            result =
                self.compute_topology(flags.contains(TriMeshFlags::DELETE_BAD_TOPOLOGY_TRIANGLES));
        }

        #[cfg(feature = "std")]
        if difference.intersects(TriMeshFlags::CONNECTED_COMPONENTS) {
            self.compute_connected_components();
        }

        #[cfg(feature = "dim3")]
        if difference.intersects(TriMeshFlags::ORIENTED | TriMeshFlags::FIX_INTERNAL_EDGES) {
            self.compute_pseudo_normals();
        }

        if prev_indices_len != self.indices.len() {
            self.rebuild_qbvh();
        }

        self.flags = flags;
        result
    }

    // TODO: support a crate like get_size2 (will require support on nalgebra too)?
    /// An approximation of the memory usage (in bytes) for this struct plus
    /// the memory it allocates dynamically.
    pub fn total_memory_size(&self) -> usize {
        size_of::<Self>() + self.heap_memory_size()
    }

    /// An approximation of the memory dynamically-allocated by this struct.
    pub fn heap_memory_size(&self) -> usize {
        // NOTE: if a new field is added to `Self`, adjust this function result.
        let Self {
            qbvh,
            vertices,
            indices,
            topology,
            connected_components,
            flags: _,
            #[cfg(feature = "dim3")]
            pseudo_normals,
        } = self;
        let sz_qbvh = qbvh.heap_memory_size();
        let sz_vertices = vertices.capacity() * size_of::<Point<Real>>();
        let sz_indices = indices.capacity() * size_of::<[u32; 3]>();
        #[cfg(feature = "dim3")]
        let sz_pseudo_normals = pseudo_normals
            .as_ref()
            .map(|pn| {
                pn.vertices_pseudo_normal.capacity() * size_of::<Vector<Real>>()
                    + pn.edges_pseudo_normal.capacity() * size_of::<[Vector<Real>; 3]>()
            })
            .unwrap_or(0);
        #[cfg(feature = "dim2")]
        let sz_pseudo_normals = 0;
        let sz_topology = topology
            .as_ref()
            .map(|t| {
                t.vertices.capacity() * size_of::<TopoVertex>()
                    + t.faces.capacity() * size_of::<TopoFace>()
                    + t.half_edges.capacity() * size_of::<TopoHalfEdge>()
            })
            .unwrap_or(0);
        let sz_connected_components = connected_components
            .as_ref()
            .map(|c| {
                c.face_colors.capacity() * size_of::<u32>()
                    + c.grouped_faces.capacity() * size_of::<f32>()
                    + c.ranges.capacity() * size_of::<usize>()
            })
            .unwrap_or(0);

        sz_qbvh
            + sz_vertices
            + sz_indices
            + sz_pseudo_normals
            + sz_topology
            + sz_connected_components
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
            pseudo_normals.edges_pseudo_normal.iter_mut().for_each(|n| {
                n[0] = transform * n[0];
                n[1] = transform * n[1];
                n[2] = transform * n[2];
            });
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
            pn.edges_pseudo_normal.iter_mut().for_each(|n| {
                n[0].component_mul_assign(scale);
                n[1].component_mul_assign(scale);
                n[2].component_mul_assign(scale);

                let _ = n[0].try_normalize_mut(0.0);
                let _ = n[1].try_normalize_mut(0.0);
                let _ = n[2].try_normalize_mut(0.0);
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

        let vertices = core::mem::take(&mut self.vertices);
        let indices = core::mem::take(&mut self.indices);
        *self = TriMesh::with_flags(vertices, indices, self.flags).unwrap();
    }

    /// Create a `TriMesh` from a set of points assumed to describe a counter-clockwise non-convex polygon.
    ///
    /// This operation may fail if the input polygon is invalid, e.g. it is non-simple or has zero surface area.
    #[cfg(feature = "dim2")]
    pub fn from_polygon(vertices: Vec<Point<Real>>) -> Option<Self> {
        triangulate_ear_clipping(&vertices).map(|indices| Self::new(vertices, indices).unwrap())
    }

    /// A flat view of the index buffer of this mesh.
    pub fn flat_indices(&self) -> &[u32] {
        unsafe {
            let len = self.indices.len() * 3;
            let data = self.indices.as_ptr() as *const u32;
            core::slice::from_raw_parts(data, len)
        }
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

        // NOTE: the Qbvh, and connected components are not changed by this operation.
        //       The pseudo-normals just have to be flipped.
        //       The topology must be recomputed.

        #[cfg(feature = "dim3")]
        if let Some(pseudo_normals) = &mut self.pseudo_normals {
            for n in &mut pseudo_normals.vertices_pseudo_normal {
                *n = -*n;
            }

            for n in pseudo_normals.edges_pseudo_normal.iter_mut() {
                n[0] = -n[0];
                n[1] = -n[1];
                n[2] = -n[2];
            }
        }

        if self.flags.contains(TriMeshFlags::HALF_EDGE_TOPOLOGY) {
            // TODO: this could be done more efficiently.
            let _ = self.compute_topology(false);
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
        let mut triangle_set = HashSet::default();

        fn resolve_coord_id(
            coord: &Point<Real>,
            vtx_to_id: &mut HashMap<HashablePartialEq<Point<Real>>, u32>,
            new_vertices: &mut Vec<Point<Real>>,
        ) -> u32 {
            let key = HashablePartialEq::new(*coord);
            let id = match vtx_to_id.entry(key) {
                Entry::Occupied(entry) => entry.into_mut(),
                Entry::Vacant(entry) => entry.insert(new_vertices.len() as u32),
            };

            if *id == new_vertices.len() as u32 {
                new_vertices.push(*coord);
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
            let _ = self.compute_topology(false);
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
    /// If the trimesh is correctly oriented, but is manifold everywhere except at its boundaries,
    /// then the computed pseudo-normals will provide correct point-containment test results except
    /// for points closest to the boundary of the mesh.
    ///
    /// It may be useful to call `self.remove_duplicate_vertices()` before this method, in order to fix the
    /// index buffer if some of the vertices of this trimesh are duplicated.
    fn compute_pseudo_normals(&mut self) {
        let mut vertices_pseudo_normal = vec![Vector::zeros(); self.vertices().len()];
        let mut edges_pseudo_normal = HashMap::default();
        let mut edges_multiplicity = HashMap::default();

        for idx in self.indices() {
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

        let edges_pseudo_normal = self
            .indices()
            .iter()
            .map(|idx| {
                let e0 = SortedPair::new(idx[0], idx[1]);
                let e1 = SortedPair::new(idx[1], idx[2]);
                let e2 = SortedPair::new(idx[2], idx[0]);
                let default = Vector::zeros();
                [
                    edges_pseudo_normal.get(&e0).copied().unwrap_or(default),
                    edges_pseudo_normal.get(&e1).copied().unwrap_or(default),
                    edges_pseudo_normal.get(&e2).copied().unwrap_or(default),
                ]
            })
            .collect();

        self.pseudo_normals = Some(TriMeshPseudoNormals {
            vertices_pseudo_normal,
            edges_pseudo_normal,
        })
    }

    fn delete_bad_topology_triangles(&mut self) {
        let mut half_edge_set = HashSet::default();
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
    fn compute_topology(&mut self, delete_bad_triangles: bool) -> Result<(), TopologyError> {
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
                    return Err(TopologyError::BadAdjacentTrianglesOrientation {
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

        Ok(())
    }

    // NOTE: this is private because that calculation is controlled by TriMeshFlags::CONNECTED_COMPONENTS
    // TODO: we should remove the CONNECTED_COMPONENTS flags and just have this be a free function.
    // TODO: this should be no_std compatible once ena is or once we have an alternative for it.
    #[cfg(feature = "std")]
    fn compute_connected_components(&mut self) {
        use ena::unify::{InPlaceUnificationTable, UnifyKey};

        #[derive(Copy, Clone, Debug, Hash, PartialEq, Eq)]
        struct IntKey(u32);

        impl UnifyKey for IntKey {
            type Value = ();
            fn index(&self) -> u32 {
                self.0
            }
            fn from_index(u: u32) -> IntKey {
                IntKey(u)
            }
            fn tag() -> &'static str {
                "IntKey"
            }
        }

        let mut ufind: InPlaceUnificationTable<IntKey> = InPlaceUnificationTable::new();
        let mut face_colors = vec![u32::MAX; self.indices.len()];
        let mut ranges = vec![0];
        let mut vertex_to_range = vec![u32::MAX; self.vertices.len()];
        let mut grouped_faces = vec![u32::MAX; self.indices.len()];
        let mut vertex_to_key = vec![IntKey(u32::MAX); self.vertices.len()];

        let mut vertex_key = |id: u32, ufind: &mut InPlaceUnificationTable<IntKey>| {
            if vertex_to_key[id as usize].0 == u32::MAX {
                let new_key = ufind.new_key(());
                vertex_to_key[id as usize] = new_key;
                new_key
            } else {
                vertex_to_key[id as usize]
            }
        };

        for idx in self.indices() {
            let keys = idx.map(|i| vertex_key(i, &mut ufind));
            ufind.union(keys[0], keys[1]);
            ufind.union(keys[1], keys[2]);
            ufind.union(keys[2], keys[0]);
        }

        for (idx, face_color) in self.indices().iter().zip(face_colors.iter_mut()) {
            debug_assert_eq!(
                ufind.find(vertex_to_key[idx[0] as usize]),
                ufind.find(vertex_to_key[idx[1] as usize])
            );
            debug_assert_eq!(
                ufind.find(vertex_to_key[idx[0] as usize]),
                ufind.find(vertex_to_key[idx[2] as usize])
            );

            let group_index = ufind.find(vertex_to_key[idx[0] as usize]).0 as usize;

            if vertex_to_range[group_index] == u32::MAX {
                // Additional range
                ranges.push(0);
                vertex_to_range[group_index] = ranges.len() as u32 - 1;
            }

            let range_id = vertex_to_range[group_index];
            ranges[range_id as usize] += 1;
            // NOTE: the range_id points to the range upper bound. The face color is the range lower bound.
            *face_color = range_id - 1;
        }

        // Cumulated sum on range indices, to find the first index faces need to be inserted into
        // for each range.
        for i in 1..ranges.len() {
            ranges[i] += ranges[i - 1];
        }

        debug_assert_eq!(*ranges.last().unwrap(), self.indices().len());

        // Group faces.
        let mut insertion_in_range_index = ranges.clone();
        for (face_id, face_color) in face_colors.iter().enumerate() {
            let insertion_index = &mut insertion_in_range_index[*face_color as usize];
            grouped_faces[*insertion_index] = face_id as u32;
            *insertion_index += 1;
        }

        self.connected_components = Some(TriMeshConnectedComponents {
            face_colors,
            grouped_faces,
            ranges,
        })
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

    #[cfg(feature = "dim3")]
    /// Gets the normal of the triangle represented by `feature`.
    pub fn feature_normal(&self, feature: FeatureId) -> Option<Unit<Vector<Real>>> {
        match feature {
            FeatureId::Face(i) => self
                .triangle(i % self.num_triangles() as u32)
                .feature_normal(FeatureId::Face(0)),
            _ => None,
        }
    }
}

impl TriMesh {
    /// The flags of this triangle mesh.
    pub fn flags(&self) -> TriMeshFlags {
        self.flags
    }

    /// Compute the axis-aligned bounding box of this triangle mesh.
    pub fn aabb(&self, pos: &Isometry<Real>) -> Aabb {
        self.qbvh.root_aabb().transform_by(pos)
    }

    /// Gets the local axis-aligned bounding box of this triangle mesh.
    pub fn local_aabb(&self) -> &Aabb {
        self.qbvh.root_aabb()
    }

    /// The acceleration structure used by this triangle-mesh.
    pub fn qbvh(&self) -> &Qbvh<u32> {
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

    /// Get the `i`-th triangle of this mesh.
    pub fn triangle(&self, i: u32) -> Triangle {
        let idx = self.indices[i as usize];
        Triangle::new(
            self.vertices[idx[0] as usize],
            self.vertices[idx[1] as usize],
            self.vertices[idx[2] as usize],
        )
    }

    /// Returns the pseudo-normals of one of this mesh’s triangles, if it was computed.
    ///
    /// This returns `None` if the pseudo-normals of this triangle were not computed.
    /// To have its pseudo-normals computed, be sure to set the [`TriMeshFlags`] so that
    /// they contain the [`TriMeshFlags::FIX_INTERNAL_EDGES`] flag.
    #[cfg(feature = "dim3")]
    pub fn triangle_normal_constraints(&self, i: u32) -> Option<TrianglePseudoNormals> {
        if self.flags.contains(TriMeshFlags::FIX_INTERNAL_EDGES) {
            let triangle = self.triangle(i);
            let pseudo_normals = self.pseudo_normals.as_ref()?;
            let edges_pseudo_normals = pseudo_normals.edges_pseudo_normal[i as usize];

            // TODO: could the pseudo-normal be pre-normalized instead of having to renormalize
            //       every time we need them?
            Some(TrianglePseudoNormals {
                face: triangle.normal()?,
                edges: [
                    Unit::try_new(edges_pseudo_normals[0], 1.0e-6)?,
                    Unit::try_new(edges_pseudo_normals[1], 1.0e-6)?,
                    Unit::try_new(edges_pseudo_normals[2], 1.0e-6)?,
                ],
            })
        } else {
            None
        }
    }

    #[cfg(feature = "dim2")]
    #[doc(hidden)]
    pub fn triangle_normal_constraints(&self, _i: u32) -> Option<TrianglePseudoNormals> {
        None
    }

    /// The vertex buffer of this mesh.
    pub fn vertices(&self) -> &[Point<Real>] {
        &self.vertices
    }

    /// The index buffer of this mesh.
    pub fn indices(&self) -> &[[u32; 3]] {
        &self.indices
    }

    /// Returns the topology information of this trimesh, if it has been computed.
    pub fn topology(&self) -> Option<&TriMeshTopology> {
        self.topology.as_ref()
    }

    /// Returns the connected-component information of this trimesh, if it has been computed.
    pub fn connected_components(&self) -> Option<&TriMeshConnectedComponents> {
        self.connected_components.as_ref()
    }

    /// Returns the connected-component of this mesh.
    ///
    /// The connected-components are returned as a set of `TriMesh` build with the given `flags`.
    pub fn connected_component_meshes(
        &self,
        flags: TriMeshFlags,
    ) -> Option<Vec<Result<TriMesh, TriMeshBuilderError>>> {
        self.connected_components()
            .map(|cc| cc.to_meshes(self, flags))
    }

    /// The pseudo-normals of this triangle mesh, if they have been computed.
    #[cfg(feature = "dim3")]
    pub fn pseudo_normals(&self) -> Option<&TriMeshPseudoNormals> {
        self.pseudo_normals.as_ref()
    }

    /// The pseudo-normals of this triangle mesh, if they have been computed **and** this mesh was
    /// marked as [`TriMeshFlags::ORIENTED`].
    #[cfg(feature = "dim3")]
    pub fn pseudo_normals_if_oriented(&self) -> Option<&TriMeshPseudoNormals> {
        if self.flags.intersects(TriMeshFlags::ORIENTED) {
            self.pseudo_normals.as_ref()
        } else {
            None
        }
    }
}

#[cfg(feature = "dim3")]
impl From<crate::shape::HeightField> for TriMesh {
    fn from(heightfield: crate::shape::HeightField) -> Self {
        let (vtx, idx) = heightfield.to_trimesh();
        TriMesh::new(vtx, idx).unwrap()
    }
}

#[cfg(feature = "dim3")]
impl From<Cuboid> for TriMesh {
    fn from(cuboid: Cuboid) -> Self {
        let (vtx, idx) = cuboid.to_trimesh();
        TriMesh::new(vtx, idx).unwrap()
    }
}

impl SimdCompositeShape for TriMesh {
    fn map_part_at(
        &self,
        i: u32,
        f: &mut dyn FnMut(Option<&Isometry<Real>>, &dyn Shape, Option<&dyn NormalConstraints>),
    ) {
        let tri = self.triangle(i);
        let normals = self.triangle_normal_constraints(i);
        f(
            None,
            &tri,
            normals.as_ref().map(|n| n as &dyn NormalConstraints),
        )
    }

    fn qbvh(&self) -> &Qbvh<u32> {
        &self.qbvh
    }
}

impl TypedSimdCompositeShape for TriMesh {
    type PartShape = Triangle;
    type PartNormalConstraints = TrianglePseudoNormals;
    type PartId = u32;

    #[inline(always)]
    fn map_typed_part_at(
        &self,
        i: u32,
        mut f: impl FnMut(
            Option<&Isometry<Real>>,
            &Self::PartShape,
            Option<&Self::PartNormalConstraints>,
        ),
    ) {
        let tri = self.triangle(i);
        let pseudo_normals = self.triangle_normal_constraints(i);
        f(None, &tri, pseudo_normals.as_ref())
    }

    #[inline(always)]
    fn map_untyped_part_at(
        &self,
        i: u32,
        mut f: impl FnMut(Option<&Isometry<Real>>, &dyn Shape, Option<&dyn NormalConstraints>),
    ) {
        let tri = self.triangle(i);
        let pseudo_normals = self.triangle_normal_constraints(i);
        f(
            None,
            &tri,
            pseudo_normals.as_ref().map(|n| n as &dyn NormalConstraints),
        )
    }

    fn typed_qbvh(&self) -> &Qbvh<u32> {
        &self.qbvh
    }
}

#[cfg(test)]
mod test {
    use crate::math::{Real, Vector};
    use crate::shape::{Cuboid, TriMesh, TriMeshFlags};

    #[test]
    fn trimesh_error_empty_indices() {
        assert!(
            TriMesh::with_flags(vec![], vec![], TriMeshFlags::empty()).is_err(),
            "A triangle mesh with no triangles is invalid."
        );
    }

    #[test]
    fn connected_components() {
        let (vtx, idx) = Cuboid::new(Vector::repeat(0.5)).to_trimesh();

        // Push 10 copy of the mesh, each time pushed with an offset.
        let mut mesh = TriMesh::new(vtx.clone(), idx.clone()).unwrap();

        for i in 1..10 {
            let cc_vtx = vtx
                .iter()
                .map(|pt| pt + Vector::repeat(2.0 * i as Real))
                .collect();

            let to_append = TriMesh::new(cc_vtx, idx.clone()).unwrap();
            mesh.append(&to_append);
        }

        mesh.set_flags(TriMeshFlags::CONNECTED_COMPONENTS).unwrap();
        let connected_components = mesh.connected_components().unwrap();
        assert_eq!(connected_components.num_connected_components(), 10);

        let cc_meshes = connected_components.to_meshes(&mesh, TriMeshFlags::empty());

        for cc in cc_meshes {
            let cc = cc.unwrap();
            assert_eq!(cc.vertices.len(), vtx.len());
            assert_eq!(cc.indices.len(), idx.len());
        }
    }
}
