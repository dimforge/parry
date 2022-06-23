//! Spatial partitioning tools.

pub use self::qbvh::{
    IndexedData, NodeIndex, QBVHDataGenerator, QbvhNonOverlappingDataSplitter, QBVH,
};
pub use self::visitor::{
    SimdBestFirstVisitStatus, SimdBestFirstVisitor, SimdSimultaneousVisitStatus,
    SimdSimultaneousVisitor, SimdVisitStatus, SimdVisitor,
};

/// A quaternary bounding-volume-hierarchy.
#[deprecated(note = "Renamed to QBVH")]
pub type SimdQbvh<T> = QBVH<T>;

mod qbvh;
mod visitor;
