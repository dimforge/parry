use na::Unit;

use crate::math::{Isometry, Point, Real, Vector};
use crate::shape::SupportMap;

/// A support mapping that is a single point.
pub struct ConstantPoint(pub Point);

impl SupportMap for ConstantPoint {
    #[inline]
    fn support_point(&self, m: &Isometry, _: &Vector) -> Point {
        m * self.0
    }

    #[inline]
    fn support_point_toward(&self, m: &Isometry, _: &Unit<Vector>) -> Point {
        m * self.0
    }

    #[inline]
    fn local_support_point(&self, _: &Vector) -> Point {
        self.0
    }

    #[inline]
    fn local_support_point_toward(&self, _: &Unit<Vector>) -> Point {
        self.0
    }
}

/// A support mapping that is the point at (0.0, 0.0, 0.0).
pub struct ConstantOrigin;

impl SupportMap for ConstantOrigin {
    #[inline]
    fn support_point(&self, m: &Isometry, _: &Vector) -> Point {
        m.translation.vector.into()
    }

    #[inline]
    fn support_point_toward(&self, m: &Isometry, _: &Unit<Vector>) -> Point {
        m.translation.vector.into()
    }

    #[inline]
    fn local_support_point(&self, _: &Vector) -> Point {
        Point::origin()
    }

    #[inline]
    fn local_support_point_toward(&self, _: &Unit<Vector>) -> Point {
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

impl<S: ?Sized + SupportMap> SupportMap for DilatedShape<'_, S> {
    #[inline]
    fn support_point(&self, m: &Isometry, dir: &Vector) -> Point {
        self.support_point_toward(m, &Unit::new_normalize(*dir))
    }

    #[inline]
    fn support_point_toward(&self, m: &Isometry, dir: &Unit<Vector>) -> Point {
        self.shape.support_point_toward(m, dir) + **dir * self.radius
    }

    #[inline]
    fn local_support_point(&self, dir: &Vector) -> Point {
        self.local_support_point_toward(&Unit::new_normalize(*dir))
    }

    #[inline]
    fn local_support_point_toward(&self, dir: &Unit<Vector>) -> Point {
        self.shape.local_support_point_toward(dir) + **dir * self.radius
    }
}
