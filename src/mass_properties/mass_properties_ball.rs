use crate::mass_properties::MassProperties;
#[cfg(feature = "dim3")]
use crate::math::Vector;
use crate::math::{Point, PrincipalAngularInertia, Real};
use na::RealField;

impl MassProperties {
    pub(crate) fn ball_volume_unit_angular_inertia(
        radius: Real,
    ) -> (Real, PrincipalAngularInertia<Real>) {
        #[cfg(feature = "dim2")]
        {
            let volume = Real::pi() * radius * radius;
            let i = radius * radius / 2.0;
            (volume, i)
        }
        #[cfg(feature = "dim3")]
        {
            let volume = Real::pi() * radius * radius * radius * 4.0 / 3.0;
            let i = radius * radius * 2.0 / 5.0;

            (volume, Vector::repeat(i))
        }
    }

    /// Computes the mass properties of a ball.
    pub fn from_ball(density: Real, radius: Real) -> Self {
        let (vol, unit_i) = Self::ball_volume_unit_angular_inertia(radius);
        let mass = vol * density;
        Self::new(Point::origin(), mass, unit_i * mass)
    }
}
