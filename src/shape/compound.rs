//!
//! Shape composed from the union of primitives.
//!

use crate::bounding_volume::{Aabb, BoundingSphere, BoundingVolume};
use crate::math::{Isometry, Real};
use crate::partitioning::Qbvh;
#[cfg(feature = "dim2")]
use crate::shape::{ConvexPolygon, TriMesh, Triangle};
use crate::shape::{Shape, SharedShape, SimdCompositeShape, TypedSimdCompositeShape};
#[cfg(feature = "dim2")]
use crate::transformation::hertel_mehlhorn;
use crate::utils::DefaultStorage;

/// A compound shape with an aabb bounding volume.
///
/// A compound shape is a shape composed of the union of several simpler shape. This is
/// the main way of creating a concave shape from convex parts. Each parts can have its own
/// delta transformation to shift or rotate it with regard to the other shapes.
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[derive(Clone, Debug)]
pub struct Compound {
    shapes: Vec<(Isometry<Real>, SharedShape)>,
    qbvh: Qbvh<u32>,
    aabbs: Vec<Aabb>,
    aabb: Aabb,
}

impl Compound {
    /// Builds a new compound shape.
    ///
    /// Panics if the input vector is empty, of if some of the provided shapes
    /// are also composite shapes (nested composite shapes are not allowed).
    pub fn new(shapes: Vec<(Isometry<Real>, SharedShape)>) -> Compound {
        assert!(
            !shapes.is_empty(),
            "A compound shape must contain at least one shape."
        );
        let mut aabbs = Vec::new();
        let mut leaves = Vec::new();
        let mut aabb = Aabb::new_invalid();

        for (i, (delta, shape)) in shapes.iter().enumerate() {
            let bv = shape.compute_aabb(delta);

            aabb.merge(&bv);
            aabbs.push(bv);
            leaves.push((i as u32, bv));

            if shape.as_composite_shape().is_some() {
                panic!("Nested composite shapes are not allowed.");
            }
        }

        let mut qbvh = Qbvh::new();
        // NOTE: we apply no dilation factor because we won't
        // update this tree dynamically.
        qbvh.clear_and_rebuild(leaves.into_iter(), 0.0);

        Compound {
            shapes,
            qbvh,
            aabbs,
            aabb,
        }
    }

    #[cfg(feature = "dim2")]
    /// Create a compound shape from the `TriMesh`. This involves merging adjacent triangles into convex
    /// polygons using the Hertel-Mehlhorn algorithm.
    ///
    /// Can fail and return `None` if any of the created shapes has close to zero or zero surface area.
    pub fn decompose_trimesh(trimesh: &TriMesh) -> Option<Self> {
        let polygons = hertel_mehlhorn(trimesh.vertices(), trimesh.indices());
        let shapes: Option<Vec<_>> = polygons
            .into_iter()
            .map(|points| {
                match points.len() {
                    3 => {
                        let triangle = Triangle::new(points[0], points[1], points[2]);
                        Some(SharedShape::new(triangle))
                    }
                    _ => ConvexPolygon::from_convex_polyline(points).map(SharedShape::new),
                }
                .map(|shape| (Isometry::identity(), shape))
            })
            .collect();
        Some(Self::new(shapes?))
    }
}

impl Compound {
    /// The shapes of this compound shape.
    #[inline]
    pub fn shapes(&self) -> &[(Isometry<Real>, SharedShape)] {
        &self.shapes[..]
    }

    /// The [`Aabb`] of this compound in its local-space.
    #[inline]
    pub fn local_aabb(&self) -> &Aabb {
        &self.aabb
    }

    /// The bounding-sphere of this compound in its local-space.
    #[inline]
    pub fn local_bounding_sphere(&self) -> BoundingSphere {
        self.aabb.bounding_sphere()
    }

    /// The shapes Aabbs.
    #[inline]
    pub fn aabbs(&self) -> &[Aabb] {
        &self.aabbs[..]
    }

    /// The acceleration structure used by this compound shape.
    #[inline]
    pub fn qbvh(&self) -> &Qbvh<u32> {
        &self.qbvh
    }
}

impl SimdCompositeShape for Compound {
    #[inline]
    fn map_part_at(&self, shape_id: u32, f: &mut dyn FnMut(Option<&Isometry<Real>>, &dyn Shape)) {
        if let Some(shape) = self.shapes.get(shape_id as usize) {
            f(Some(&shape.0), &*shape.1)
        }
    }

    #[inline]
    fn qbvh(&self) -> &Qbvh<u32> {
        &self.qbvh
    }
}

impl TypedSimdCompositeShape for Compound {
    type PartShape = dyn Shape;
    type PartId = u32;
    type QbvhStorage = DefaultStorage;

    #[inline(always)]
    fn map_typed_part_at(
        &self,
        i: u32,
        mut f: impl FnMut(Option<&Isometry<Real>>, &Self::PartShape),
    ) {
        if let Some((part_pos, part)) = self.shapes.get(i as usize) {
            f(Some(part_pos), &**part)
        }
    }

    #[inline(always)]
    fn map_untyped_part_at(
        &self,
        i: u32,
        mut f: impl FnMut(Option<&Isometry<Real>>, &Self::PartShape),
    ) {
        if let Some((part_pos, part)) = self.shapes.get(i as usize) {
            f(Some(part_pos), &**part)
        }
    }

    #[inline]
    fn typed_qbvh(&self) -> &Qbvh<u32> {
        &self.qbvh
    }
}
