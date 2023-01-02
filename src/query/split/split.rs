/// The result of a plane-splitting operation.
pub enum SplitResult<T> {
    /// The split operation yield two results: one lying on the negative half-space of the plane
    /// and the second lying on the positive half-space of the plane.
    Pair(T, T),
    /// The shape being split is fully contained in the negative half-space of the plane.
    Negative,
    /// The shape being split is fully contained in the positive half-space of the plane.
    Positive,
}

/// The result of a plane-intersection operation.
pub enum IntersectResult<T> {
    /// The intersect operation yielded a result, lying in the plane
    Intersect(T),
    /// The shape being intersected is fully contained in the negative half-space of the plane.
    Negative,
    /// The shape being intersected is fully contained in the positive half-space of the plane.
    Positive,
}
