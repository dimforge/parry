use crate::mass_properties::MassProperties;
use crate::math::{Point, Real};
#[cfg(feature = "dim3")]
use crate::shape::Capsule;

impl MassProperties {
    /// Computes the mass properties of a capsule.
    pub fn from_capsule(density: Real, a: Point<Real>, b: Point<Real>, radius: Real) -> Self {
        let half_height = (b - a).norm() / 2.0;
        let (cyl_vol, cyl_unit_i) = Self::cylinder_y_volume_unit_inertia(half_height, radius);
        let (ball_vol, ball_unit_i) = Self::ball_volume_unit_angular_inertia(radius);
        let cap_vol = cyl_vol + ball_vol;
        let cap_mass = cap_vol * density;
        let mut cap_i = (cyl_unit_i * cyl_vol + ball_unit_i * ball_vol) * density;
        let local_com = na::center(&a, &b);

        #[cfg(feature = "dim2")]
        {
            let h = half_height * 2.0;
            let extra = (h * h * 0.25 + h * radius * 3.0 / 8.0) * ball_vol * density;
            cap_i += extra;
            Self::new(local_com, cap_mass, cap_i)
        }

        #[cfg(feature = "dim3")]
        {
            let h = half_height * 2.0;
            let extra = (h * h * 0.25 + h * radius * 3.0 / 8.0) * ball_vol * density;
            cap_i.x += extra;
            cap_i.z += extra;
            let local_frame = Capsule::new(a, b, radius).rotation_wrt_y();
            Self::with_principal_inertia_frame(local_com, cap_mass, cap_i, local_frame)
        }
    }
}
