use crate::mass_properties::MassProperties;
use crate::math::{Point, PrincipalAngularInertia, Real, Vector};

impl MassProperties {
    pub(crate) fn cuboid_volume_unit_inertia(
        half_extents: Vector<Real>,
    ) -> (Real, PrincipalAngularInertia<Real>) {
        #[cfg(feature = "dim2")]
        {
            let volume = half_extents.x * half_extents.y * 4.0;
            let ix = (half_extents.x * half_extents.x) / 3.0;
            let iy = (half_extents.y * half_extents.y) / 3.0;

            (volume, ix + iy)
        }

        #[cfg(feature = "dim3")]
        {
            let volume = half_extents.x * half_extents.y * half_extents.z * 8.0;
            let ix = (half_extents.x * half_extents.x) / 3.0;
            let iy = (half_extents.y * half_extents.y) / 3.0;
            let iz = (half_extents.z * half_extents.z) / 3.0;

            (volume, Vector::new(iy + iz, ix + iz, ix + iy))
        }
    }

    /// Computes the mass properties of a cuboid.
    pub fn from_cuboid(density: Real, half_extents: Vector<Real>) -> Self {
        let (vol, unit_i) = Self::cuboid_volume_unit_inertia(half_extents);
        let mass = vol * density;
        Self::new(Point::origin(), mass, unit_i * mass)
    }
}
