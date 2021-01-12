//! Transformation, simplification and decomposition of meshes.

#[cfg(feature = "dim3")]
pub(crate) use self::convex_hull2::convex_hull2_idx;
#[cfg(feature = "dim2")]
pub use self::convex_hull2::{convex_hull2 as convex_hull, convex_hull2_idx as convex_hull_idx};
#[cfg(feature = "dim3")]
pub use self::convex_hull3::convex_hull3 as convex_hull;

mod convex_hull2;
#[cfg(feature = "dim3")]
mod convex_hull3;
pub(crate) mod convex_hull_utils;

#[cfg(feature = "dim3")]
pub mod vhacd;

#[cfg(feature = "dim2")]
mod to_polyline;
#[cfg(feature = "dim3")]
mod to_trimesh;
mod utils;
