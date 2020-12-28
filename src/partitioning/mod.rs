//! Spatial partitioning tools.

pub use self::visitor::{
    SimdBestFirstVisitStatus, SimdBestFirstVisitor, SimdSimultaneousVisitor, SimdVisitStatus,
    SimdVisitor,
};
pub use self::wquadtree::{IndexedData, SimdQuadTree};

mod visitor;
mod wquadtree;
