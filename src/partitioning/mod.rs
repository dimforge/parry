//! Spatial partitioning tools.

#[cfg(feature = "std")]
pub use self::qbvh::{CenterDataSplitter, QBVHDataGenerator, QBVHNonOverlappingDataSplitter};
pub use self::qbvh::{
    GenericQBVH, IndexedData, NodeIndex, QBVHNode, QBVHProxy, QBVHStorage, SimdNodeIndex, QBVH,
};
#[cfg(feature = "parallel")]
pub use self::visitor::{ParallelSimdSimultaneousVisitor, ParallelSimdVisitor};
pub use self::visitor::{
    SimdBestFirstVisitStatus, SimdBestFirstVisitor, SimdSimultaneousVisitStatus,
    SimdSimultaneousVisitor, SimdVisitStatus, SimdVisitor,
};

/// A quaternary bounding-volume-hierarchy.
#[deprecated(note = "Renamed to QBVH")]
pub type SimdQBVH<T> = QBVH<T>;

mod qbvh;
mod visitor;
