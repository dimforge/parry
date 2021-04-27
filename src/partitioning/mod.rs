//! Spatial partitioning tools.

pub use self::visitor::{
    SimdBestFirstVisitStatus, SimdBestFirstVisitor, SimdSimultaneousVisitor, SimdVisitStatus,
    SimdVisitor,
};
pub use self::wquadtree::{IndexedData, QBVHDataGenerator, QBVH};

#[deprecated(note = "Renamed to QBVH")]
pub type SimdQuadTree<T> = QBVH<T>;

mod visitor;
mod wquadtree;
