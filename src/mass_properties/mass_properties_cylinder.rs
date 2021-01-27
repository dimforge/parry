use crate::mass_properties::MassProperties;
use crate::math::{PrincipalAngularInertia, Real, Vector};
#[cfg(feature = "dim3")]
use {
    crate::math::{Point, Rotation},
    na::RealField,
};

impl MassProperties {
    pub(crate) fn cylinder_y_volume_unit_inertia(
        half_height: Real,
        radius: Real,
    ) -> (Real, PrincipalAngularInertia<Real>) {
        #[cfg(feature = "dim2")]
        {
            Self::cuboid_volume_unit_inertia(Vector::new(radius, half_height))
        }

        #[cfg(feature = "dim3")]
        {
            let volume = half_height * radius * radius * Real::pi() * 2.0;
            let sq_radius = radius * radius;
            let sq_height = half_height * half_height * 4.0;
            let off_principal = (sq_radius * 3.0 + sq_height) / 12.0;

            let inertia = Vector::new(off_principal, sq_radius / 2.0, off_principal);
            (volume, inertia)
        }
    }

    /// Computes the mass properties of a cylinder.
    #[cfg(feature = "dim3")]
    pub fn from_cylinder(density: Real, half_height: Real, radius: Real) -> Self {
        let (cyl_vol, cyl_unit_i) = Self::cylinder_y_volume_unit_inertia(half_height, radius);
        let cyl_mass = cyl_vol * density;

        Self::with_principal_inertia_frame(
            Point::origin(),
            cyl_mass,
            cyl_unit_i * cyl_mass,
            Rotation::identity(),
        )
    }
}
