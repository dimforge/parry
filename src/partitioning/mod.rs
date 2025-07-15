//! Spatial partitioning tools.

#[cfg(feature = "alloc")]
pub use self::bvh::{Bvh, BvhBuildStrategy, BvhLeafCost, BvhNode, BvhWorkspace, TraversalAction};

#[cfg(feature = "alloc")]
mod bvh;
