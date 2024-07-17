//! Various unsorted geometrical and logical operators.

pub use self::ccw_face_normal::ccw_face_normal;
pub use self::center::center;
#[cfg(feature = "std")]
pub use self::deterministic_state::DeterministicState;

#[cfg(feature = "dim3")]
#[cfg(feature = "std")]
pub use self::cleanup::remove_unused_points;
pub(crate) use self::inv::inv;
pub use self::isometry_ops::{IsometryOps, IsometryOpt};
pub use self::median::median;
pub use self::point_cloud_support_point::{
    point_cloud_support_point, point_cloud_support_point_id,
};
pub use self::point_in_poly2d::{point_in_convex_poly2d, point_in_poly2d};
pub use self::sdp_matrix::{SdpMatrix2, SdpMatrix3};

pub use self::as_bytes::AsBytes;
pub(crate) use self::consts::*;
pub use self::cov::{center_cov, cov};
pub use self::hashable_partial_eq::HashablePartialEq;
#[cfg(feature = "std")]
pub use self::interval::{find_root_intervals, find_root_intervals_to, Interval, IntervalFunction};
pub use self::obb::obb;
pub use self::segments_intersection::{segments_intersection2d, SegmentsIntersection};
#[cfg(feature = "dim3")]
pub(crate) use self::sort::sort2;
pub(crate) use self::sort::sort3;
pub use self::sorted_pair::SortedPair;
#[cfg(all(feature = "dim3", feature = "std"))]
pub(crate) use self::spade::sanitize_spade_point;
pub(crate) use self::weighted_value::WeightedValue;
pub(crate) use self::wops::{simd_swap, WBasis, WCross, WSign};

mod as_bytes;
mod ccw_face_normal;
mod center;
#[cfg(feature = "dim3")]
#[cfg(feature = "std")]
mod cleanup;
mod consts;
mod cov;
#[cfg(feature = "std")]
mod deterministic_state;
mod hashable_partial_eq;
#[cfg(feature = "std")]
pub mod hashmap;
#[cfg(feature = "std")]
mod interval;
mod inv;
mod isometry_ops;
mod median;
mod obb;
mod point_cloud_support_point;
mod point_in_poly2d;
#[cfg(feature = "dim2")]
pub mod point_in_triangle;
mod sdp_matrix;
mod segments_intersection;
mod sort;
mod sorted_pair;
#[cfg(all(feature = "dim3", feature = "std"))]
mod spade;
mod weighted_value;
mod wops;
