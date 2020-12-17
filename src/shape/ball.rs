use na::Unit;

use crate::math::{Isometry, Point, Real, Vector};
use crate::shape::SupportMap;

/// A Ball shape.
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[derive(PartialEq, Debug, Copy, Clone)]
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
