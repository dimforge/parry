//! Support mapping based HalfSpace shape.
use crate::math::{Real, Vector};
use na::Unit;

/// A half-space delimited by an infinite plane.
#[derive(PartialEq, Debug, Clone, Copy)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[cfg_attr(
    all(not(target_os = "cuda"), feature = "cuda"),
    derive(cust::DeviceCopy)
)]
#[repr(C)]
pub struct HalfSpace {
    /// The halfspace planar boundary outward normal.
    pub normal: Unit<Vector<Real>>,
}

impl HalfSpace {
    /// Builds a new halfspace from its center and its normal.
    #[inline]
    pub fn new(normal: Unit<Vector<Real>>) -> HalfSpace {
        HalfSpace { normal }
    }
}
