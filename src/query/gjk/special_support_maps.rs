use na::Unit;

use crate::math::{Isometry, Point, Real, Vector};
use crate::shape::SupportMap;

/// A support mapping that is a single point.
pub struct ConstantPoint(pub Point<Real>);

impl SupportMap for ConstantPoint {
    #[inline]
    fn support_point(&self, m: &Isometry<Real>, _: &Vector<Real>) -> Point<Real> {
        m * self.0
    }

    #[inline]
    fn support_point_toward(&self, m: &Isometry<Real>, _: &Unit<Vector<Real>>) -> Point<Real> {
        m * self.0
    }

    #[inline]
    fn local_support_point(&self, _: &Vector<Real>) -> Point<Real> {
        self.0
    }

    #[inline]
    fn local_support_point_toward(&self, _: &Unit<Vector<Real>>) -> Point<Real> {
        self.0
    }
}

/// A support mapping that is the point at (0.0, 0.0, 0.0).
pub struct ConstantOrigin;

impl SupportMap for ConstantOrigin {
    #[inline]
    fn support_point(&self, m: &Isometry<Real>, _: &Vector<Real>) -> Point<Real> {
        m.translation.vector.into()
    }

    #[inline]
    fn support_point_toward(&self, m: &Isometry<Real>, _: &Unit<Vector<Real>>) -> Point<Real> {
        m.translation.vector.into()
    }

    #[inline]
    fn local_support_point(&self, _: &Vector<Real>) -> Point<Real> {
        Point::origin()
    }

    #[inline]
    fn local_support_point_toward(&self, _: &Unit<Vector<Real>>) -> Point<Real> {
        Point::origin()
    }
}

/// The Minkowski sum of a shape and a ball.
pub struct DilatedShape<'a, S: ?Sized + SupportMap> {
    /// The shape involved in the Minkowski sum.
    pub shape: &'a S,
    /// The radius of the ball involved in the Minkoski sum.
    pub radius: Real,
}

impl<'a, S: ?Sized + SupportMap> SupportMap for DilatedShape<'a, S> {
    #[inline]
    fn support_point(&self, m: &Isometry<Real>, dir: &Vector<Real>) -> Point<Real> {
        self.support_point_toward(m, &Unit::new_normalize(*dir))
    }

    #[inline]
    fn support_point_toward(&self, m: &Isometry<Real>, dir: &Unit<Vector<Real>>) -> Point<Real> {
        self.shape.support_point_toward(m, dir) + **dir * self.radius
    }

    #[inline]
    fn local_support_point(&self, dir: &Vector<Real>) -> Point<Real> {
        self.local_support_point_toward(&Unit::new_normalize(*dir))
    }

    #[inline]
    fn local_support_point_toward(&self, dir: &Unit<Vector<Real>>) -> Point<Real> {
        self.shape.local_support_point_toward(dir) + **dir * self.radius
    }
}
