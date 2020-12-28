use crate::math::{Isometry, Real};
use crate::partitioning::SimdQuadTree;
use crate::shape::Shape;

/// Trait implemented by shapes composed of multiple simpler shapes.
///
/// A composite shape is composed of several shapes. For example, this can
/// be a convex decomposition of a concave shape; or a triangle-mesh.
pub trait SimdCompositeShape {
    /// The number of sub-shape in this composite shape.
    fn nparts(&self) -> usize;

    /// Applies a function to one sub-shape of this composite shape.
    fn map_part_at(&self, shape_id: u32, f: &mut dyn FnMut(Option<&Isometry<Real>>, &dyn Shape));

    /// Gets the acceleration structure of the composite shape.
    fn quadtree(&self) -> &SimdQuadTree<u32>;
}

pub trait TypedSimdCompositeShape: SimdCompositeShape {
    type PartShape: ?Sized + Shape;

    fn map_typed_part_at(
        &self,
        shape_id: u32,
        f: impl FnMut(Option<&Isometry<Real>>, &Self::PartShape),
    );
}

impl TypedSimdCompositeShape for dyn SimdCompositeShape {
    type PartShape = dyn Shape;
    fn map_typed_part_at(
        &self,
        shape_id: u32,
        mut f: impl FnMut(Option<&Isometry<Real>>, &Self::PartShape),
    ) {
        self.map_part_at(shape_id, &mut f)
    }
}
