use crate::mass_properties::MassProperties;
use crate::math::{Point, Real, real};
use crate::shape::Triangle;

impl MassProperties {
    /// Computes the mass properties of a triangle.
    pub fn from_triangle(
        density: Real,
        a: &Point<Real>,
        b: &Point<Real>,
        c: &Point<Real>,
    ) -> MassProperties {
        let triangle = Triangle::new(*a, *b, *c);
        let area = triangle.area();
        let com = triangle.center();

        if area == real!(0.0) {
            return MassProperties::new(com, real!(0.0), real!(0.0));
        }

        let ipart = triangle.unit_angular_inertia();

        Self::new(com, area * density, ipart * area * density)
    }
}
