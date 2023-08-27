#[cfg(feature = "std")]
use either::Either;
use na::Unit;

use crate::math::{Isometry, Point, Real, Vector};
use crate::shape::SupportMap;

/// A Ball shape.
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[cfg_attr(feature = "bytemuck", derive(bytemuck::Pod, bytemuck::Zeroable))]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Deserialize, rkyv::Serialize),
    archive(check_bytes)
)]
#[cfg_attr(feature = "cuda", derive(cust_core::DeviceCopy))]
#[derive(PartialEq, Debug, Copy, Clone)]
#[repr(C)]
pub struct Ball {
    /// The radius of the ball.
    pub radius: Real,
}

impl Ball {
    /// Creates a new ball with the given radius.
    #[inline]
    pub fn new(radius: Real) -> Ball {
        Ball { radius }
    }

    /// Computes a scaled version of this ball.
    ///
    /// If the scaling factor is non-uniform, then it can’t be represented as
    /// ball. Instead, a convex polygon approximation (with `nsubdivs`
    /// subdivisions) is returned. Returns `None` if that approximation had degenerate
    /// normals (for example if the scaling factor along one axis is zero).
    #[cfg(all(feature = "dim2", feature = "std"))]
    #[inline]
    pub fn scaled(
        self,
        scale: &Vector<Real>,
        nsubdivs: u32,
    ) -> Option<Either<Self, super::ConvexPolygon>> {
        if scale.x != scale.y {
            // The scaled shape isn’t a ball.
            let mut vtx = self.to_polyline(nsubdivs);
            vtx.iter_mut()
                .for_each(|pt| pt.coords = pt.coords.component_mul(scale));
            Some(Either::Right(super::ConvexPolygon::from_convex_polyline(
                vtx,
            )?))
        } else {
            let uniform_scale = scale.x;
            Some(Either::Left(Self::new(self.radius * uniform_scale.abs())))
        }
    }

    /// Computes a scaled version of this ball.
    ///
    /// If the scaling factor is non-uniform, then it can’t be represented as
    /// ball. Instead, a convex polygon approximation (with `nsubdivs`
    /// subdivisions) is returned. Returns `None` if that approximation had degenerate
    /// normals (for example if the scaling factor along one axis is zero).
    #[cfg(all(feature = "dim3", feature = "std"))]
    #[inline]
    pub fn scaled(
        self,
        scale: &Vector<Real>,
        nsubdivs: u32,
    ) -> Option<Either<Self, super::ConvexPolyhedron>> {
        if scale.x != scale.y || scale.x != scale.z || scale.y != scale.z {
            // The scaled shape isn’t a ball.
            let (mut vtx, idx) = self.to_trimesh(nsubdivs, nsubdivs);
            vtx.iter_mut()
                .for_each(|pt| pt.coords = pt.coords.component_mul(scale));
            Some(Either::Right(super::ConvexPolyhedron::from_convex_mesh(
                vtx, &idx,
            )?))
        } else {
            let uniform_scale = scale.x;
            Some(Either::Left(Self::new(self.radius * uniform_scale.abs())))
        }
    }
}

impl SupportMap for Ball {
    #[inline]
    fn support_point(&self, m: &Isometry<Real>, dir: &Vector<Real>) -> Point<Real> {
        self.support_point_toward(m, &Unit::new_normalize(*dir))
    }

    #[inline]
    fn support_point_toward(&self, m: &Isometry<Real>, dir: &Unit<Vector<Real>>) -> Point<Real> {
        Point::from(m.translation.vector) + **dir * self.radius
    }

    #[inline]
    fn local_support_point(&self, dir: &Vector<Real>) -> Point<Real> {
        self.local_support_point_toward(&Unit::new_normalize(*dir))
    }

    #[inline]
    fn local_support_point_toward(&self, dir: &Unit<Vector<Real>>) -> Point<Real> {
        Point::from(**dir * self.radius)
    }
}
