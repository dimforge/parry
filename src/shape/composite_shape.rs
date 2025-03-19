use crate::math::Isometry;
use crate::partitioning::{IndexedData, Qbvh};
use crate::query::details::NormalConstraints;
use crate::shape::Shape;

/// Trait implemented by shapes composed of multiple simpler shapes.
///
/// A composite shape is composed of several shapes. For example, this can
/// be a convex decomposition of a concave shape; or a triangle-mesh.
#[cfg(feature = "std")]
pub trait SimdCompositeShape {
    /// Applies a function to one sub-shape of this composite shape.
    fn map_part_at(
        &self,
        shape_id: u32,
        f: &mut dyn FnMut(Option<&Isometry>, &dyn Shape, Option<&dyn NormalConstraints>),
    );

    /// Gets the acceleration structure of the composite shape.
    fn qbvh(&self) -> &Qbvh<u32>;
}

#[cfg(feature = "std")]
pub trait TypedSimdCompositeShape {
    type PartShape: ?Sized + Shape;
    type PartNormalConstraints: ?Sized + NormalConstraints;
    type PartId: IndexedData;

    fn map_typed_part_at(
        &self,
        shape_id: Self::PartId,
        f: impl FnMut(Option<&Isometry>, &Self::PartShape, Option<&Self::PartNormalConstraints>),
    );

    // TODO: we need this method because the compiler won't want
    // to cast `&Self::PartShape` to `&dyn Shape` because it complains
    // that `PartShape` is not `Sized`.
    fn map_untyped_part_at(
        &self,
        shape_id: Self::PartId,
        f: impl FnMut(Option<&Isometry>, &dyn Shape, Option<&dyn NormalConstraints>),
    );

    fn typed_qbvh(&self) -> &Qbvh<Self::PartId>;
}

#[cfg(feature = "std")]
impl TypedSimdCompositeShape for dyn SimdCompositeShape + '_ {
    type PartShape = dyn Shape;
    type PartNormalConstraints = dyn NormalConstraints;
    type PartId = u32;

    fn map_typed_part_at(
        &self,
        shape_id: u32,
        mut f: impl FnMut(Option<&Isometry>, &Self::PartShape, Option<&Self::PartNormalConstraints>),
    ) {
        self.map_part_at(shape_id, &mut f)
    }

    fn map_untyped_part_at(
        &self,
        shape_id: u32,
        mut f: impl FnMut(Option<&Isometry>, &dyn Shape, Option<&dyn NormalConstraints>),
    ) {
        self.map_part_at(shape_id, &mut f)
    }

    fn typed_qbvh(&self) -> &Qbvh<Self::PartId> {
        self.qbvh()
    }
}
