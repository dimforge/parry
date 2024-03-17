use crate::math::*;
use crate::shape::SupportMap;

/// A support mapping that is a single point.
pub struct ConstantPoint(pub Point);

impl SupportMap for ConstantPoint {
    #[inline]
    fn support_point(&self, m: &Isometry, _: &Vector) -> Point {
        m.transform_point(&self.0)
    }

    #[inline]
    fn support_point_toward(&self, m: &Isometry, _: &UnitVector) -> Point {
        m.transform_point(&self.0)
    }

    #[inline]
    fn local_support_point(&self, _: &Vector) -> Point {
        self.0
    }

    #[inline]
    fn local_support_point_toward(&self, _: &UnitVector) -> Point {
        self.0
    }
}

/// A support mapping that is the point at (0.0, 0.0, 0.0).
pub struct ConstantOrigin;

impl SupportMap for ConstantOrigin {
    #[inline]
    fn support_point(&self, m: &Isometry, _: &Vector) -> Point {
        m.translation.into()
    }

    #[inline]
    fn support_point_toward(&self, m: &Isometry, _: &UnitVector) -> Point {
        m.translation.into()
    }

    #[inline]
    fn local_support_point(&self, _: &Vector) -> Point {
        Point::origin()
    }

    #[inline]
    fn local_support_point_toward(&self, _: &UnitVector) -> Point {
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
    fn support_point(&self, m: &Isometry, dir: &Vector) -> Point {
        self.support_point_toward(m, &UnitVector::new_normalize(*dir))
    }

    #[inline]
    fn support_point_toward(&self, m: &Isometry, dir: &UnitVector) -> Point {
        self.shape.support_point_toward(m, dir) + dir.into_inner() * self.radius
    }

    #[inline]
    fn local_support_point(&self, dir: &Vector) -> Point {
        self.local_support_point_toward(&UnitVector::new_normalize(*dir))
    }

    #[inline]
    fn local_support_point_toward(&self, dir: &UnitVector) -> Point {
        self.shape.local_support_point_toward(dir) + dir.into_inner() * self.radius
    }
}
