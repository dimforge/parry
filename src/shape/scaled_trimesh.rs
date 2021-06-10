use std::sync::Arc;

use crate::bounding_volume::AABB;
use crate::math::{Isometry, Point, Real, Vector};
use crate::partitioning::QBVH;
use crate::shape::composite_shape::SimdCompositeShape;
use crate::shape::TriMesh;
use crate::shape::{Shape, Triangle, TypedSimdCompositeShape};

#[derive(Clone)]
#[cfg_attr(feature = "serde-serialize", derive(Serialize, Deserialize))]
/// A scaled [`TriMesh`]
pub struct ScaledTriMesh {
    /// The underlying triangle mesh
    trimesh: Arc<TriMesh>,
    /// Scaling factors for each dimension
    // This could easily be expanded into an arbitrary transform, if a use case arises.
    scaling_factors: Vector<Real>,
    quadtree: QBVH<u32>,
}

impl ScaledTriMesh {
    /// Creates a triangle mesh by scaling `trimesh` along each axis
    pub fn new(trimesh: Arc<TriMesh>, scaling_factors: Vector<Real>) -> Self {
        // Future work: Would it be more efficient to scale trimesh.quadtree rather than building a
        // new one from scratch?
        let data = trimesh.triangles().enumerate().map(|(i, tri)| {
            let aabb = scale_tri(&tri, &scaling_factors).local_aabb();
            (i as u32, aabb)
        });
        let mut quadtree = QBVH::new();
        quadtree.clear_and_rebuild(data, 0.0);
        Self {
            trimesh,
            scaling_factors,
            quadtree,
        }
    }

    /// The underlying unscaled trimesh
    pub fn trimesh(&self) -> &Arc<TriMesh> {
        &self.trimesh
    }

    /// Scaling factors used to derive this shape from the underlying trimesh
    pub fn scaling_factors(&self) -> Vector<Real> {
        self.scaling_factors
    }

    /// Compute the axis-aligned bounding box of this triangle mesh.
    pub fn aabb(&self, pos: &Isometry<Real>) -> AABB {
        self.quadtree.root_aabb().transform_by(pos)
    }

    /// Gets the local axis-aligned bounding box of this triangle mesh.
    pub fn local_aabb(&self) -> &AABB {
        self.quadtree.root_aabb()
    }

    /// The acceleration structure used by this triangle-mesh.
    pub fn quadtree(&self) -> &QBVH<u32> {
        &self.quadtree
    }

    /// An iterator through all the scaled triangles of this mesh.
    pub fn triangles(&self) -> impl ExactSizeIterator<Item = Triangle> + '_ {
        self.trimesh
            .triangles()
            .map(move |tri| scale_tri(&tri, &self.scaling_factors))
    }

    /// Get the `i`-th scaled triangle of this mesh.
    pub fn triangle(&self, i: u32) -> Triangle {
        scale_tri(&self.trimesh.triangle(i), &self.scaling_factors)
    }

    /// The vertex buffer of this mesh.
    pub fn vertices(&self) -> impl ExactSizeIterator<Item = Point<Real>> + '_ {
        self.trimesh
            .vertices()
            .iter()
            .map(move |v| v.coords.component_mul(&self.scaling_factors).into())
    }

    /// The index buffer of this mesh.
    pub fn indices(&self) -> &[[u32; 3]] {
        self.trimesh.indices()
    }
}

fn scale_tri(tri: &Triangle, factors: &Vector<Real>) -> Triangle {
    Triangle {
        a: tri.a.coords.component_mul(factors).into(),
        b: tri.b.coords.component_mul(factors).into(),
        c: tri.c.coords.component_mul(factors).into(),
    }
}

impl SimdCompositeShape for ScaledTriMesh {
    fn map_part_at(&self, i: u32, f: &mut dyn FnMut(Option<&Isometry<Real>>, &dyn Shape)) {
        let tri = self.triangle(i);
        f(None, &tri)
    }

    fn quadtree(&self) -> &QBVH<u32> {
        &self.quadtree
    }
}

impl TypedSimdCompositeShape for ScaledTriMesh {
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

    fn typed_quadtree(&self) -> &QBVH<u32> {
        &self.quadtree
    }
}
