use crate::math::{Isometry, Real};
use crate::partitioning::{GenericQbvh, IndexedData, Qbvh, QbvhStorage};
use crate::shape::Shape;
use crate::utils::DefaultStorage;

/// Trait implemented by shapes composed of multiple simpler shapes.
///
/// A composite shape is composed of several shapes. For example, this can
/// be a convex decomposition of a concave shape; or a triangle-mesh.
#[cfg(feature = "std")]
pub trait SimdCompositeShape {
    /// Applies a function to one sub-shape of this composite shape.
    fn map_part_at(&self, shape_id: u32, f: &mut dyn FnMut(Option<&Isometry<Real>>, &dyn Shape));

    /// Gets the acceleration structure of the composite shape.
    fn qbvh(&self) -> &Qbvh<u32>;
}

pub trait TypedSimdCompositeShape {
    type PartShape: ?Sized + Shape;
    type PartId: IndexedData;
    type QbvhStorage: QbvhStorage<Self::PartId>;

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

    fn typed_qbvh(&self) -> &GenericQbvh<Self::PartId, Self::QbvhStorage>;
}

#[cfg(feature = "std")]
impl<'a> TypedSimdCompositeShape for dyn SimdCompositeShape + 'a {
    type PartShape = dyn Shape;
    type PartId = u32;
    type QbvhStorage = DefaultStorage;

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

    fn typed_qbvh(&self) -> &GenericQbvh<Self::PartId, Self::QbvhStorage> {
        self.qbvh()
    }
}
