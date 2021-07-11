//! Spatial partitioning tools.

pub use self::qbvh::{IndexedData, QBVHDataGenerator, QBVH};
pub use self::visitor::{
    SimdBestFirstVisitStatus, SimdBestFirstVisitor, SimdSimultaneousVisitor, SimdVisitStatus,
    SimdVisitor,
};

/// A quaternary bounding-volume-hierarchy.
#[deprecated(note = "Renamed to QBVH")]
pub type SimdQbvh<T> = QBVH<T>;

mod qbvh;
mod visitor;
