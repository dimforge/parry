//! Spatial partitioning tools.

#[cfg(feature = "alloc")]
pub use self::bvh::{Bvh, BvhBuildStrategy, BvhLeafCost, BvhNode, BvhWorkspace};

#[cfg(feature = "alloc")]
mod bvh;
