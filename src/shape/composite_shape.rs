use crate::partitioning::WQuadtree;
use crate::shape::Shape;

/// Trait implemented by shapes composed of multiple simpler shapes.
///
/// A composite shape is composed of several shapes. Typically, it is a convex decomposition of
/// a concave shape.
pub trait SimdCompositeShape {
    /// The number of sub-shape in this composite shape.
    fn nparts(&self) -> usize;

    /// Applies a transformation matrix and a function to each sub-shape of this composite
    /// shape.
    fn map_part_at(&self, _: u32, _: &mut dyn FnMut(&dyn Shape));

    /// Gets the acceleration structure of the composite shape.
    fn quadtree(&self) -> &WQuadtree<u32>;
}
