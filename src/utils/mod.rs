//! Various unsorted geometrical and logical operators.

pub use self::ccw_face_normal::ccw_face_normal;
pub use self::center::center;
pub use self::deterministic_state::DeterministicState;

#[cfg(feature = "dim3")]
pub use self::cleanup::remove_unused_points;
pub(crate) use self::inv::inv;
pub use self::isometry_ops::{IsometryOps, IsometryOpt};
pub use self::median::median;
pub use self::point_cloud_support_point::{
    point_cloud_support_point, point_cloud_support_point_id,
};
pub use self::point_in_poly2d::point_in_poly2d;
pub use self::sdp_matrix::{SdpMatrix2, SdpMatrix3};

pub use self::as_bytes::AsBytes;
pub(crate) use self::consts::*;
pub use self::cov::cov;
pub use self::hashable_partial_eq::HashablePartialEq;
pub use self::interval::{find_root_intervals, find_root_intervals_to, Interval, IntervalFunction};
#[cfg(feature = "dim3")]
pub(crate) use self::sort::sort2;
pub(crate) use self::sort::sort3;
pub use self::sorted_pair::SortedPair;
pub(crate) use self::weighted_value::WeightedValue;
pub(crate) use self::wops::{simd_swap, WBasis, WCross, WSign};

mod as_bytes;
mod ccw_face_normal;
mod center;
#[cfg(feature = "dim3")]
mod cleanup;
mod consts;
mod cov;
mod deterministic_state;
mod hashable_partial_eq;
pub mod hashmap;
mod interval;
mod inv;
mod isometry_ops;
mod median;
mod point_cloud_support_point;
mod point_in_poly2d;
mod ref_with_cost;
mod sdp_matrix;
mod sort;
mod sorted_pair;
mod weighted_value;
mod wops;
