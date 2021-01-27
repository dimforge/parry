//! Transformation, simplification and decomposition of meshes.

#[cfg(feature = "dim3")]
pub(crate) use self::convex_hull2::convex_hull2_idx;
#[cfg(feature = "dim2")]
pub use self::convex_hull2::{convex_hull2 as convex_hull, convex_hull2_idx as convex_hull_idx};
#[cfg(feature = "dim3")]
pub use self::convex_hull3::{check_convex_hull, convex_hull};

mod convex_hull2;
#[cfg(feature = "dim3")]
mod convex_hull3;
pub(crate) mod convex_hull_utils;

/// Approximate convex decomposition using the VHACD algorithm.
pub mod vhacd;
/// Voxelization of a 2D polyline or 3D triangle mesh.
pub mod voxelization;

#[cfg(feature = "dim2")]
mod to_polyline;
#[cfg(feature = "dim3")]
mod to_trimesh;
mod utils;
