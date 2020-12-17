use crate::math::{Isometry, Point, Real, Vector};
use crate::shape::SupportMap;
use na::Unit;
use std::ops::Sub;

/// A point of a Configuration-Space Obstacle.
///
/// A Configuration-Space Obstacle (CSO) is the result of the
/// Minkowski Difference of two solids. In other words, each of its
/// points correspond to the difference of two point, each belonging
/// to a different solid.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct CSOPoint {
    /// The point on the CSO. This is equal to `self.orig1 - self.orig2`, unless this CSOPoint
    /// has been translated with self.translate.
    pub point: Point<Real>,
    /// The original point on the first shape used to compute `self.point`.
    pub orig1: Point<Real>,
    /// The original point on the second shape used to compute `self.point`.
    pub orig2: Point<Real>,
}

impl CSOPoint {
    /// Initializes a CSO point with `orig1 - orig2`.
    pub fn new(orig1: Point<Real>, orig2: Point<Real>) -> Self {
        let point = Point::from(orig1 - orig2);
        Self::new_with_point(point, orig1, orig2)
    }

    /// Initializes a CSO point with all information provided.
    ///
    /// It is assumed, but not checked, that `point == orig1 - orig2`.
    pub fn new_with_point(point: Point<Real>, orig1: Point<Real>, orig2: Point<Real>) -> Self {
        CSOPoint {
            point,
            orig1,
            orig2,
        }
    }

    /// Initializes a CSO point where both original points are equal.
    pub fn single_point(point: Point<Real>) -> Self {
        Self::new_with_point(point, point, Point::origin())
    }

    /// CSO point where all components are set to zero.
    pub fn origin() -> Self {
        CSOPoint::new(Point::origin(), Point::origin())
    }

    /// Computes the support point of the CSO of `g1` and `g2` toward the unit direction `dir`.
    pub fn from_shapes_toward<G1: ?Sized, G2: ?Sized>(
        pos12: &Isometry<Real>,
        g1: &G1,
        g2: &G2,
        dir: &Unit<Vector<Real>>,
    ) -> Self
    where
        G1: SupportMap,
        G2: SupportMap,
    {
        let sp1 = g1.local_support_point_toward(dir);
        let sp2 = g2.support_point_toward(pos12, &-*dir);

        CSOPoint::new(sp1, sp2)
    }

    /// Computes the support point of the CSO of `g1` and `g2` toward the direction `dir`.
    pub fn from_shapes<G1: ?Sized, G2: ?Sized>(
        pos12: &Isometry<Real>,
        g1: &G1,
        g2: &G2,
        dir: &Vector<Real>,
    ) -> Self
    where
        G1: SupportMap,
        G2: SupportMap,
    {
        let sp1 = g1.local_support_point(dir);
        let sp2 = g2.support_point(pos12, &-*dir);

        CSOPoint::new(sp1, sp2)
    }

    /// Translate the CSO point.
    pub fn translate(&self, dir: &Vector<Real>) -> Self {
        CSOPoint::new_with_point(self.point + dir, self.orig1, self.orig2)
    }

    /// Translate in-place the CSO point.
    pub fn translate_mut(&mut self, dir: &Vector<Real>) {
        self.point += dir;
    }
}

impl Sub<CSOPoint> for CSOPoint {
    type Output = Vector<Real>;

    #[inline]
    fn sub(self, rhs: CSOPoint) -> Vector<Real> {
        self.point - rhs.point
    }
}
