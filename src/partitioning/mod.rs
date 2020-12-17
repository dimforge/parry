//! Spatial partitioning tools.

pub use self::visitor::{
    SimdBestFirstVisitStatus, SimdBestFirstVisitor, SimdSimultaneousVisitor, SimdVisitStatus,
    SimdVisitor,
};
pub use self::wquadtree::{IndexedData, WQuadtree};

mod visitor;
mod wquadtree;
