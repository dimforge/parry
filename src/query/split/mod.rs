pub use self::split::{CanonicalSplit, Split, SplitResult};

mod split;
mod split_aabb;
mod split_convex_polygon;
mod split_convex_polyhedron;
mod split_polyline;
mod split_segment;

#[cfg(all(feature = "std", feature = "dim3"))]
mod split_trimesh;
