//! Spatial partitioning tools.

pub use self::qbvh::{
    CenterDataSplitter, IndexedData, NodeIndex, QBVHDataGenerator, QBVHNode, QBVHNodeDataGenerator,
    QBVHProxy, QbvhNonOverlappingDataSplitter, SimdNodeIndex, QBVH,
};
#[cfg(feature = "parallel")]
pub use self::visitor::{ParallelSimdSimultaneousVisitor, ParallelSimdVisitor};
pub use self::visitor::{
    SimdBestFirstVisitStatus, SimdBestFirstVisitor, SimdSimultaneousVisitStatus,
    SimdSimultaneousVisitor, SimdVisitStatus, SimdVisitor,
};

/// A quaternary bounding-volume-hierarchy.
#[deprecated(note = "Renamed to QBVH")]
pub type SimdQbvh<T> = QBVH<T>;

mod qbvh;
mod visitor;
