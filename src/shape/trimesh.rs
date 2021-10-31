use crate::bounding_volume::AABB;
use crate::math::{Isometry, Point, Real};
use crate::partitioning::QBVH;
use crate::shape::composite_shape::SimdCompositeShape;
#[cfg(feature = "dim2")]
use crate::shape::{Compound, ConvexPolygon, SharedShape};
use crate::shape::{FeatureId, Shape, Triangle, TypedSimdCompositeShape};
#[cfg(feature = "dim2")]
use crate::transformation::ear_clipping::triangulate_ear_clipping;
#[cfg(feature = "dim2")]
use crate::transformation::hertel_mehlhorn::hertel_mehlhorn;
use crate::utils::hashmap::{Entry, HashMap};
use crate::utils::HashablePartialEq;
#[cfg(feature = "dim3")]
use {
    crate::math::Vector,
    crate::shape::{Cuboid, HeightField},
    crate::utils::SortedPair,
};

#[derive(Clone)]
#[cfg_attr(feature = "serde-serialize", derive(Serialize, Deserialize))]
#[cfg(feature = "dim3")]
pub(crate) struct PseudoNormals {
    pub vertices_pseudo_normal: Vec<Vector<Real>>,
    pub edges_pseudo_normal: HashMap<SortedPair<u32>, Vector<Real>>,
}

#[cfg(feature = "dim3")]
impl Default for PseudoNormals {
    fn default() -> Self {
        Self {
            vertices_pseudo_normal: vec![],
            edges_pseudo_normal: Default::default(),
        }
    }
}

#[derive(Clone)]
#[cfg_attr(feature = "serde-serialize", derive(Serialize, Deserialize))]
/// A triangle mesh.
pub struct TriMesh {
    qbvh: QBVH<u32>,
    vertices: Vec<Point<Real>>,
    indices: Vec<[u32; 3]>,
    #[cfg(feature = "dim3")]
    pub(crate) pseudo_normals: PseudoNormals,
}

impl TriMesh {
    /// Creates a new triangle mesh from a vertex buffer and an index buffer.
    pub fn new(vertices: Vec<Point<Real>>, indices: Vec<[u32; 3]>) -> Self {
        assert!(
            indices.len() > 0,
            "A triangle mesh must contain at least one triangle."
        );

        let data = indices.iter().enumerate().map(|(i, idx)| {
            let aabb = Triangle::new(
                vertices[idx[0] as usize],
                vertices[idx[1] as usize],
                vertices[idx[2] as usize],
            )
            .local_aabb();
            (i as u32, aabb)
        });

        let mut qbvh = QBVH::new();
        // NOTE: we apply no dilation factor because we won't
        // update this tree dynamically.
        qbvh.clear_and_rebuild(data, 0.0);

        Self {
            qbvh,
            vertices,
            indices,
            #[cfg(feature = "dim3")]
            pseudo_normals: Default::default(),
        }
    }

    /// Create a `TriMesh` from a set of points assumed to describe a counter-clockwise non-convex polygon.
    ///
    /// This operation may fail if the input polygon is invalid, e.g. it is non-simple or has zero surface area.
    #[cfg(feature = "dim2")]
    pub fn from_polygon(vertices: Vec<Point<Real>>) -> Option<Self> {
        triangulate_ear_clipping(&vertices).map(|indices| Self::new(vertices, indices))
    }

    #[cfg(feature = "dim2")]
    /// Create a compound shape from the `TriMesh`. This involves merging adjacent triangles into convex
    /// polygons using the Hertel-Mehlhorn algorithm.
    ///
    /// Can fail and return `None` if any of the created shapes has close to zero or zero surface area.
    pub fn to_compound(&self) -> Option<Compound> {
        let indices = hertel_mehlhorn(self.vertices(), self.indices());
        let shapes: Option<Vec<_>> = indices
            .into_iter()
            .map(|poly_indices| {
                match poly_indices.as_slice() {
                    [i0, i1, i2] => {
                        let triangle = Triangle::new(
                            self.vertices()[*i0 as usize],
                            self.vertices()[*i1 as usize],
                            self.vertices()[*i2 as usize],
                        );
                        Some(SharedShape::new(triangle))
                    }
                    i_polygon => {
                        let vertices = i_polygon
                            .iter()
                            .map(|i| self.vertices()[*i as usize])
                            .collect();
                        ConvexPolygon::from_convex_polyline(vertices).map(SharedShape::new)
                    }
                }
                .map(|shape| (Isometry::identity(), shape))
            })
            .collect();
        Some(Compound::new(shapes?))
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
    pub fn triangles(&self) -> impl Iterator<Item = Triangle> + '_ {
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

    /// Remove all duplicate vertices and adjust the index buffer accordingly.
    ///
    /// This is typically used to recover a vertex buffer from which we can deduce
    /// adjacency information. between triangles by observing how the vertices are
    /// shared by triangles based on the index buffer.
    pub fn recover_topology(&mut self) {
        let mut vtx_to_id = HashMap::default();
        let mut new_vertices = Vec::with_capacity(self.vertices.len());
        let mut new_indices = Vec::with_capacity(self.indices.len());

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

            new_indices.push([va, vb, vc]);
        }

        new_vertices.shrink_to_fit();

        self.vertices = new_vertices;
        self.indices = new_indices;

        // Vertices and indices changed: the pseudo-normals are no longer valid.
        #[cfg(feature = "dim3")]
        if !self.pseudo_normals.vertices_pseudo_normal.is_empty() {
            self.compute_pseudo_normals();
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
    /// It may be useful to call `self.recover_topology()` before this method, in order to fix the
    /// index buffer if some of the vertices of this trimesh are duplicated.
    pub fn compute_pseudo_normals(&mut self) {
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

        self.pseudo_normals.vertices_pseudo_normal = vertices_pseudo_normal;
        self.pseudo_normals.edges_pseudo_normal = edges_pseudo_normal;
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
impl From<HeightField> for TriMesh {
    fn from(heightfield: HeightField) -> Self {
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
