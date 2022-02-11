use crate::math::{Isometry, Real, UnitVector};

pub enum SplitResult<T> {
    Pair(T, T),
    Negative,
    Positive,
}

pub trait CanonicalSplit: Sized {
    /// Splits this shape along the given canonical axis.
    ///
    /// This will split the shape by a plane with a normal with it’s `axis`-th component set to 1.
    /// The splittin plan is shifted wrt. the origin by the `bias` (i.e. it passes through the point
    /// equal to `normal * bias`).
    ///
    /// # Result
    /// Returns the result of the split. The first sahpe returned is the piece lying on the negative
    /// half-space delimited by the splitting plane. The second sahpe returned is the piece lying on the
    /// positive half-space delimited by the splitting plane.
    fn canonical_split(&self, axis: usize, bias: Real, eps: Real) -> SplitResult<Self>;
}

pub trait Split: Sized + CanonicalSplit {
    fn split(
        &self,
        position: &Isometry<Real>,
        axis: &UnitVector<Real>,
        bias: Real,
        epsilon: Real,
    ) -> SplitResult<Self> {
        let additional_bias = position.translation.vector.dot(axis);
        self.local_split(
            &position.inverse_transform_unit_vector(axis),
            bias + additional_bias,
            epsilon,
        )
    }

    fn local_split(
        &self,
        local_axis: &UnitVector<Real>,
        bias: Real,
        epsilon: Real,
    ) -> SplitResult<Self>;
}

// TODO: impl CanonicalSplit for TriMesh, Polyline, ConvexPolygon, ConvexPolyhedron, Cuboid.
// TODO: impl Split for TriMesh, Polyline, ConvexPolygon, ConvexPolyhedron.
