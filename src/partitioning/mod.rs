//! Spatial partitioning tools.

#[cfg(feature = "alloc")]
pub use self::qbvh::{
    CenterDataSplitter, IndexedData, NodeIndex, Qbvh, QbvhDataGenerator, QbvhNode,
    QbvhNonOverlappingDataSplitter, QbvhProxy, QbvhUpdateWorkspace, SimdNodeIndex,
};
#[cfg(feature = "parallel")]
pub use self::visitor::{ParallelSimdSimultaneousVisitor, ParallelSimdVisitor};
pub use self::visitor::{
    SimdBestFirstVisitStatus, SimdBestFirstVisitor, SimdSimultaneousVisitStatus,
    SimdSimultaneousVisitor, SimdVisitStatus, SimdVisitor, SimdVisitorWithContext,
};

/// A quaternary bounding-volume-hierarchy.
#[deprecated(note = "Renamed to Qbvh")]
#[cfg(feature = "alloc")]
pub type SimdQbvh<T> = Qbvh<T>;

#[cfg(feature = "alloc")]
mod qbvh;
mod visitor;
