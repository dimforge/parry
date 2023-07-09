//! Support mapping based HalfSpace shape.
use crate::math::{Real, Vector};
use na::Unit;

#[cfg(feature = "rkyv")]
use rkyv::{bytecheck, CheckBytes};

/// A half-space delimited by an infinite plane.
#[derive(PartialEq, Debug, Clone, Copy)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Deserialize, rkyv::Serialize, CheckBytes),
    archive(as = "Self")
)]
#[cfg_attr(feature = "cuda", derive(cust_core::DeviceCopy))]
#[repr(C)]
pub struct HalfSpace {
    /// The halfspace planar boundary's outward normal.
    pub normal: Unit<Vector<Real>>,
}

impl HalfSpace {
    /// Builds a new halfspace from its center and its normal.
    #[inline]
    pub fn new(normal: Unit<Vector<Real>>) -> HalfSpace {
        HalfSpace { normal }
    }

    /// Computes a scaled version of this half-space.
    ///
    /// Returns `None` if `self.normal` scaled by `scale` is zero (the scaled half-space
    /// degenerates to a single point).
    pub fn scaled(self, scale: &Vector<Real>) -> Option<Self> {
        Unit::try_new(self.normal.component_mul(scale), 0.0).map(|normal| Self { normal })
    }
}
