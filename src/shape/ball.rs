use either::Either;
use na::Unit;

use crate::math::{Isometry, Point, Real, Vector};
use crate::shape::SupportMap;

/// A Ball shape.
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[cfg_attr(feature = "cuda", derive(cust_core::DeviceCopy))]
#[derive(PartialEq, Debug, Copy, Clone)]
#[repr(C)]
pub struct Ball {
    /// The radius of the ball.
    pub radius: Real,
}

impl Ball {
    /// Creates a new ball from its radius and center.
    #[inline]
    pub fn new(radius: Real) -> Ball {
        Ball { radius }
    }

    #[cfg(feature = "dim2")]
    #[inline]
    pub fn scaled(
        self,
        scale: &Vector<Real>,
        nsubdivs: u32,
    ) -> Option<Either<Self, super::ConvexPolygon>> {
        if scale.x != scale.y {
            // The scaled shape isn’t a ball.
            let mut vtx = Self::new(0.5).to_polyline(nsubdivs);
            vtx.iter_mut()
                .for_each(|pt| pt.coords = pt.coords.component_mul(&scale));
            Some(Either::Right(super::ConvexPolygon::from_convex_polyline(
                vtx,
            )?))
        } else {
            let uniform_scale = scale.x;
            Some(Either::Left(Self::new(self.radius * uniform_scale.abs())))
        }
    }

    #[cfg(feature = "dim3")]
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
                .for_each(|pt| pt.coords = pt.coords.component_mul(&scale));
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
