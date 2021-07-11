//!
//! Shape composed from the union of primitives.
//!

use crate::bounding_volume::{BoundingSphere, BoundingVolume, AABB};
use crate::math::{Isometry, Real};
use crate::partitioning::QBVH;
use crate::shape::{Shape, SharedShape, SimdCompositeShape, TypedSimdCompositeShape};

/// A compound shape with an aabb bounding volume.
///
/// A compound shape is a shape composed of the union of several simpler shape. This is
/// the main way of creating a concave shape from convex parts. Each parts can have its own
/// delta transformation to shift or rotate it with regard to the other shapes.
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[derive(Clone)]
pub struct Compound {
    shapes: Vec<(Isometry<Real>, SharedShape)>,
    qbvh: QBVH<u32>,
    aabbs: Vec<AABB>,
    aabb: AABB,
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
        let mut aabb = AABB::new_invalid();

        for (i, &(ref delta, ref shape)) in shapes.iter().enumerate() {
            let bv = shape.compute_aabb(delta);

            aabb.merge(&bv);
            aabbs.push(bv.clone());
            leaves.push((i as u32, bv));

            if shape.as_composite_shape().is_some() {
                panic!("Nested composite shapes are not allowed.");
            }
        }

        let mut qbvh = QBVH::new();
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
}

impl Compound {
    /// The shapes of this compound shape.
    #[inline]
    pub fn shapes(&self) -> &[(Isometry<Real>, SharedShape)] {
        &self.shapes[..]
    }

    /// The AABB of this compound in its local-space.
    #[inline]
    pub fn local_aabb(&self) -> &AABB {
        &self.aabb
    }

    /// The bounding-sphere of this compound in its local-space.
    #[inline]
    pub fn local_bounding_sphere(&self) -> BoundingSphere {
        self.aabb.bounding_sphere()
    }

    /// The shapes AABBs.
    #[inline]
    pub fn aabbs(&self) -> &[AABB] {
        &self.aabbs[..]
    }

    /// The acceleration structure used by this compound shape.
    #[inline]
    pub fn qbvh(&self) -> &QBVH<u32> {
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
    fn qbvh(&self) -> &QBVH<u32> {
        &self.qbvh
    }
}

impl TypedSimdCompositeShape for Compound {
    type PartShape = dyn Shape;
    type PartId = u32;

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
    fn typed_qbvh(&self) -> &QBVH<u32> {
        &self.qbvh
    }
}
