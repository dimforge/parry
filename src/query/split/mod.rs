pub use self::split::{IntersectResult, SplitResult};

mod split;
mod split_aabb;
mod split_segment;

#[cfg(all(any(feature = "std", feature = "alloc"), feature = "dim3"))]
mod split_trimesh;
