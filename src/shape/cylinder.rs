//! Support mapping based Cylinder shape.

use crate::math::{Point, Real, Vector};
use crate::shape::SupportMap;
use na;
use num::Zero;

/// Cylinder shape with its principal axis aligned with the `y` axis.
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[derive(PartialEq, Debug, Copy, Clone)]
pub struct Cylinder {
    /// The half-height of the cylinder.
    pub half_height: Real,
    /// The radius fo the cylinder.
    pub radius: Real,
}

impl Cylinder {
    /// Creates a new cylinder.
    ///
    /// # Arguments:
    /// * `half_height` - the half length of the cylinder along the `y` axis.
    /// * `radius` - the length of the cylinder along all other axis.
    pub fn new(half_height: Real, radius: Real) -> Cylinder {
        assert!(half_height.is_sign_positive() && radius.is_sign_positive());

        Cylinder {
            half_height,
            radius,
        }
    }
}

impl SupportMap for Cylinder {
    fn local_support_point(&self, dir: &Vector<Real>) -> Point<Real> {
        let mut vres = *dir;

        vres[1] = 0.0;

        if vres.normalize_mut().is_zero() {
            vres = na::zero()
        } else {
            vres = vres * self.radius;
        }

        vres[1] = self.half_height.copysign(dir[1]);

        Point::from(vres)
    }
}
