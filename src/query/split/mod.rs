pub use self::split::{IntersectResult, SplitResult};

mod split;
mod split_aabb;
mod split_segment;

#[cfg(all(feature = "dim3", feature = "spade"))]
mod split_trimesh;
