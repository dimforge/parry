use crate::math::{Isometry, Real};
use crate::partitioning::Bvh;
use crate::query::details::NormalConstraints;
use crate::shape::Shape;

/// Trait implemented by shapes composed of multiple simpler shapes.
///
/// A composite shape is composed of several shapes. For example, this can
/// be a convex decomposition of a concave shape; or a triangle-mesh.
#[cfg(feature = "alloc")]
pub trait SimdCompositeShape {
    /// Applies a function to one sub-shape of this composite shape.
    fn map_part_at(
        &self,
        shape_id: u32,
        f: &mut dyn FnMut(Option<&Isometry<Real>>, &dyn Shape, Option<&dyn NormalConstraints>),
    );

    /// Gets the acceleration structure of the composite shape.
    fn bvh(&self) -> &Bvh;
}

#[cfg(feature = "alloc")]
pub trait TypedSimdCompositeShape {
    type PartShape: ?Sized + Shape;
    type PartNormalConstraints: ?Sized + NormalConstraints;

    // TODO: have this return a T instead of taking an FnMut?
    fn map_typed_part_at<T>(
        &self,
        shape_id: u32,
        f: impl FnMut(
            Option<&Isometry<Real>>,
            &Self::PartShape,
            Option<&Self::PartNormalConstraints>,
        ) -> T,
    ) -> Option<T>;

    // TODO: we need this method because the compiler won't want
    // to cast `&Self::PartShape` to `&dyn Shape` because it complains
    // that `PartShape` is not `Sized`.
    fn map_untyped_part_at<T>(
        &self,
        shape_id: u32,
        f: impl FnMut(Option<&Isometry<Real>>, &dyn Shape, Option<&dyn NormalConstraints>) -> T,
    ) -> Option<T>;

    fn typed_bvh(&self) -> &Bvh;
}

#[cfg(feature = "alloc")]
impl TypedSimdCompositeShape for dyn SimdCompositeShape + '_ {
    type PartShape = dyn Shape;
    type PartNormalConstraints = dyn NormalConstraints;

    fn map_typed_part_at<T>(
        &self,
        shape_id: u32,
        mut f: impl FnMut(
            Option<&Isometry<Real>>,
            &Self::PartShape,
            Option<&Self::PartNormalConstraints>,
        ) -> T,
    ) -> Option<T> {
        let mut result = None;
        self.map_part_at(shape_id, &mut |pose, part, normals| {
            result = Some(f(pose, part, normals));
        });
        result
    }

    fn map_untyped_part_at<T>(
        &self,
        shape_id: u32,
        mut f: impl FnMut(Option<&Isometry<Real>>, &dyn Shape, Option<&dyn NormalConstraints>) -> T,
    ) -> Option<T> {
        let mut result = None;
        self.map_part_at(shape_id, &mut |pose, part, normals| {
            result = Some(f(pose, part, normals));
        });
        result
    }

    fn typed_bvh(&self) -> &Bvh {
        self.bvh()
    }
}

/// A helper struct that implements scene queries on any composite shapes.
///
/// For example, the `RayCast` implementation of a composite shape can use this wrapper or
/// provide its own implementation. This is for working around the lack of specialization in
/// (stable) rust. If we did have specialization, this would just be a blanket implementation
/// of all the geometric query traits for all `S: CompositeShape`.
pub struct CompositeShapeRef<'a, S: ?Sized>(pub &'a S);
