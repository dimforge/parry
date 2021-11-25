//! Support mapping based Cone shape.

use crate::math::{Point, Real, Vector};
use crate::shape::SupportMap;
use na;
use num::Zero;

#[cfg(not(feature = "std"))]
use na::RealField; // for .copysign()

/// Cone shape with its principal axis aligned with the `y` axis.
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[cfg_attr(
    all(not(target_os = "cuda"), feature = "cuda"),
    derive(cust::DeviceCopy)
)]
#[derive(PartialEq, Debug, Copy, Clone)]
#[repr(C)]
pub struct Cone {
    /// The half-height of the cone.
    pub half_height: Real,
    /// The base radius of the cone.
    pub radius: Real,
}

impl Cone {
    /// Creates a new cone.
    ///
    /// # Arguments:
    /// * `half_height` - the half length of the cone along the `y` axis.
    /// * `radius` - the length of the cone along all other axis.
    pub fn new(half_height: Real, radius: Real) -> Cone {
        Cone {
            half_height,
            radius,
        }
    }
}

impl SupportMap for Cone {
    #[inline]
    fn local_support_point(&self, dir: &Vector<Real>) -> Point<Real> {
        let mut vres = *dir;

        vres[1] = 0.0;

        if vres.normalize_mut().is_zero() {
            vres = na::zero();
            vres[1] = self.half_height.copysign(dir[1]);
        } else {
            vres = vres * self.radius;
            vres[1] = -self.half_height;

            if dir.dot(&vres) < dir[1] * self.half_height {
                vres = na::zero();
                vres[1] = self.half_height
            }
        }

        Point::from(vres)
    }
}
