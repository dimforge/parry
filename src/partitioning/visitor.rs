use crate::math::{Real, SimdBool, SimdReal, SIMD_WIDTH};

/// The next action to be taken by a BVH traversal algorithm after having visited a node with some data.
pub enum SimdBestFirstVisitStatus<Res> {
    MaybeContinue {
        weights: SimdReal,
        mask: SimdBool,
        results: [Option<Res>; SIMD_WIDTH],
    },
    /// The traversal aborts.
    ///
    /// If a data is provided, then it is returned as the result of the traversal.
    /// If no result is provided, then the last best result found becomes the result of the traversal.
    ExitEarly(Option<Res>),
}

/// Trait implemented by cost functions used by the best-first search on a `BVT`.
pub trait SimdBestFirstVisitor<T, SimdBV> {
    /// The result of a best-first traversal.
    type Result;

    /// Compute the next action to be taken by the best-first-search after visiting a node containing the given bounding volume.
    fn visit(
        &mut self,
        best_cost_so_far: Real,
        bv: &SimdBV,
        value: Option<[Option<&T>; SIMD_WIDTH]>,
    ) -> SimdBestFirstVisitStatus<Self::Result>;
}

/// The status of the spatial partitioning structure traversal.
pub enum SimdVisitStatus {
    /// The traversal should continue on the children of the currently visited nodes for which
    /// the boolean lane is set to `1`.
    MaybeContinue(SimdBool),
    /// The traversal should exit immediately.
    ExitEarly,
}

/// Trait implemented by visitor called during the traversal of a spatial partitioning data structure.
pub trait SimdVisitor<T, SimdBV> {
    /// Execute an operation on the content of a node of the spatial partitioning structure.
    ///
    /// Returns whether the traversal should continue on the node's children, if it should not continue
    /// on those children, or if the whole traversal should be exited early.
    fn visit(&mut self, bv: &SimdBV, data: Option<[Option<&T>; SIMD_WIDTH]>) -> SimdVisitStatus;
}

/// Trait implemented by visitor called during a simultaneous spatial partitioning data structure tarversal.
pub trait SimdSimultaneousVisitor<T, SimdBV> {
    /// Execute an operation on the content of two nodes, one from each structure.
    ///
    /// Returns whether the traversal should continue on the nodes children, if it should not continue
    /// on those children, or if the whole traversal should be exited early.
    fn visit(
        &mut self,
        left_bv: &SimdBV,
        left_data: Option<[Option<&T>; SIMD_WIDTH]>,
        right_bv: &SimdBV,
        right_data: Option<[Option<&T>; SIMD_WIDTH]>,
    ) -> SimdVisitStatus;
}
