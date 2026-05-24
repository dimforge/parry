//! This module only exists for working around the fact that glam
//! doesn’t allow its `approx` feature to be enabled when targeting spirv.

use crate::math::Vector;
use approx::relative_eq;

/// Checks approximate equality between two vectors, component-wise.
#[inline]
pub fn relative_eq_vector(a: Vector, b: Vector) -> bool {
    #[cfg(feature = "dim2")]
    return relative_eq!(a.x, b.x) && relative_eq!(a.y, b.y);
    #[cfg(feature = "dim3")]
    return relative_eq!(a.x, b.x) && relative_eq!(a.y, b.y) && relative_eq!(a.z, b.z);
}
