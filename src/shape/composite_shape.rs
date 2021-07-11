use crate::math::{Isometry, Real};
use crate::partitioning::{IndexedData, QBVH};
use crate::shape::Shape;

/// Trait implemented by shapes composed of multiple simpler shapes.
///
/// A composite shape is composed of several shapes. For example, this can
/// be a convex decomposition of a concave shape; or a triangle-mesh.
pub trait SimdCompositeShape {
    /// Applies a function to one sub-shape of this composite shape.
    fn map_part_at(&self, shape_id: u32, f: &mut dyn FnMut(Option<&Isometry<Real>>, &dyn Shape));

    /// Gets the acceleration structure of the composite shape.
    fn qbvh(&self) -> &QBVH<u32>;
}

pub trait TypedSimdCompositeShape {
    type PartShape: ?Sized + Shape;
    type PartId: IndexedData;

    fn map_typed_part_at(
        &self,
        shape_id: Self::PartId,
        f: impl FnMut(Option<&Isometry<Real>>, &Self::PartShape),
    );

    // TODO: we need this method because the compiler won't want
    // to cast `&Self::PartShape` to `&dyn Shape` because it complains
    // that `PairtShape` is not `Sized`.
    fn map_untyped_part_at(
        &self,
        shape_id: Self::PartId,
        f: impl FnMut(Option<&Isometry<Real>>, &dyn Shape),
    );

    fn typed_qbvh(&self) -> &QBVH<Self::PartId>;
}

impl<'a> TypedSimdCompositeShape for dyn SimdCompositeShape + 'a {
    type PartShape = dyn Shape;
    type PartId = u32;

    fn map_typed_part_at(
        &self,
        shape_id: u32,
        mut f: impl FnMut(Option<&Isometry<Real>>, &Self::PartShape),
    ) {
        self.map_part_at(shape_id, &mut f)
    }

    fn map_untyped_part_at(
        &self,
        shape_id: u32,
        mut f: impl FnMut(Option<&Isometry<Real>>, &dyn Shape),
    ) {
        self.map_part_at(shape_id, &mut f)
    }

    fn typed_qbvh(&self) -> &QBVH<Self::PartId> {
        self.qbvh()
    }
}
