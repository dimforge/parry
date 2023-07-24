use crate::mass_properties::MassProperties;
use crate::math::{Point, PrincipalAngularInertia, Real, real, Rotation, Vector};
use na::RealField;

impl MassProperties {
    pub(crate) fn cone_y_volume_unit_inertia(
        half_height: Real,
        radius: Real,
    ) -> (Real, PrincipalAngularInertia<Real>) {
        let volume = radius * radius * Real::pi() * half_height * real!(2.0) / real!(3.0);
        let sq_radius = radius * radius;
        let sq_height = half_height * half_height * real!(4.0);
        let off_principal = sq_radius * real!(3.0) / real!(20.0) + sq_height * real!(3.0) / real!(80.0);
        let principal = sq_radius * real!(3.0) / real!(10.0);

        (volume, Vector::new(off_principal, principal, off_principal))
    }

    /// Computes the mass properties of a cone.
    pub fn from_cone(density: Real, half_height: Real, radius: Real) -> Self {
        let (cyl_vol, cyl_unit_i) = Self::cone_y_volume_unit_inertia(half_height, radius);
        let cyl_mass = cyl_vol * density;

        Self::with_principal_inertia_frame(
            Point::new(real!(0.0), -half_height / real!(2.0), real!(0.0)),
            cyl_mass,
            cyl_unit_i * cyl_mass,
            Rotation::identity(),
        )
    }
}
