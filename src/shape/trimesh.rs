use crate::bounding_volume::AABB;
use crate::math::{Isometry, Point, Real};
use crate::partitioning::SimdQuadTree;
use crate::shape::composite_shape::SimdCompositeShape;
#[cfg(feature = "dim3")]
use crate::shape::{Cuboid, HeightField};
use crate::shape::{Shape, Triangle, TypedSimdCompositeShape};

#[derive(Clone)]
#[cfg_attr(feature = "serde-serialize", derive(Serialize, Deserialize))]
/// A triangle mesh.
pub struct TriMesh {
    quadtree: SimdQuadTree<u32>,
    vertices: Vec<Point<Real>>,
    indices: Vec<[u32; 3]>,
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

        let mut quadtree = SimdQuadTree::new();
        // NOTE: we apply no dilation factor because we won't
        // update this tree dynamically.
        quadtree.clear_and_rebuild(data, 0.0);

        Self {
            quadtree,
            vertices,
            indices,
        }
    }

    /// Compute the axis-aligned bounding box of this triangle mesh.
    pub fn aabb(&self, pos: &Isometry<Real>) -> AABB {
        self.quadtree.root_aabb().transform_by(pos)
    }

    /// Gets the local axis-aligned bounding box of this triangle mesh.
    pub fn local_aabb(&self) -> &AABB {
        self.quadtree.root_aabb()
    }

    pub fn quadtree(&self) -> &SimdQuadTree<u32> {
        &self.quadtree
    }

    /// The number of triangles forming this mesh.
    pub fn num_triangles(&self) -> usize {
        self.indices.len()
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
        self.quadtree.cast_ray(&ray, max_toi, &mut intersections);
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
        self.quadtree.cast_ray(&ray, max_toi, &mut intersections);

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

    fn quadtree(&self) -> &SimdQuadTree<u32> {
        &self.quadtree
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

    fn typed_quadtree(&self) -> &SimdQuadTree<u32> {
        &self.quadtree
    }
}
