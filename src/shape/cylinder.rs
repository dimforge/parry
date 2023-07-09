//! Support mapping based Cylinder shape.

use crate::math::{Point, Real, Vector};
use crate::shape::SupportMap;
use na;
use num::Zero;

#[cfg(feature = "std")]
use either::Either;

#[cfg(not(feature = "std"))]
use na::RealField; // for .copysign()

#[cfg(feature = "rkyv")]
use rkyv::{bytecheck, CheckBytes};

/// Cylinder shape with its principal axis aligned with the `y` axis.
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[cfg_attr(feature = "bytemuck", derive(bytemuck::Pod, bytemuck::Zeroable))]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Deserialize, rkyv::Serialize, CheckBytes),
    archive(as = "Self")
)]
#[cfg_attr(feature = "cuda", derive(cust_core::DeviceCopy))]
#[derive(PartialEq, Debug, Copy, Clone)]
#[repr(C)]
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

    /// Computes a scaled version of this cylinder.
    ///
    /// If the scaling factor is non-uniform, then it can’t be represented as
    /// cylinder. Instead, a convex polyhedral approximation (with `nsubdivs`
    /// subdivisions) is returned. Returns `None` if that approximation had degenerate
    /// normals (for example if the scaling factor along one axis is zero).
    #[cfg(feature = "std")]
    #[inline]
    pub fn scaled(
        self,
        scale: &Vector<Real>,
        nsubdivs: u32,
    ) -> Option<Either<Self, super::ConvexPolyhedron>> {
        if scale.x != scale.z {
            // The scaled shape isn’t a cylinder.
            let (mut vtx, idx) = self.to_trimesh(nsubdivs);
            vtx.iter_mut()
                .for_each(|pt| pt.coords = pt.coords.component_mul(&scale));
            Some(Either::Right(super::ConvexPolyhedron::from_convex_mesh(
                vtx, &idx,
            )?))
        } else {
            Some(Either::Left(Self::new(
                self.half_height * scale.y,
                self.radius * scale.x,
            )))
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
