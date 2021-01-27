use crate::math::{Real, SimdBool, SimdReal, SIMD_WIDTH};

/// The next action to be taken by a BVH traversal algorithm after having visited a node with some data.
pub enum SimdBestFirstVisitStatus<Res> {
    /// The traversal can continue.
    MaybeContinue {
        /// The weight associated to each child of the node being traversed.
        weights: SimdReal,
        /// Each lane indicates if the corresponding child of the node being traversed
        /// should be traversed too.
        mask: SimdBool,
        /// Optional results associated to each child of the node being traversed.
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

/// The status of the simultaneous traversal of two spatial partitioning structures.
pub enum SimdSimultaneousVisitStatus {
    /// The traversal should continue on the children of the currently visited nodes for which
    /// the boolean lane is set to `1`.
    MaybeContinue([SimdBool; SIMD_WIDTH]),
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

impl<F, T, SimdBV> SimdVisitor<T, SimdBV> for F
where
    F: FnMut(&SimdBV, Option<[Option<&T>; SIMD_WIDTH]>) -> SimdVisitStatus,
{
    fn visit(&mut self, bv: &SimdBV, data: Option<[Option<&T>; SIMD_WIDTH]>) -> SimdVisitStatus {
        (self)(bv, data)
    }
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
    ) -> SimdSimultaneousVisitStatus;
}
