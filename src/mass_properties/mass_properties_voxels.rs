use crate::mass_properties::MassProperties;
use crate::math::{Point, Real};
use crate::shape::Voxels;

impl MassProperties {
    /// Computes the mass properties of a set of voxels.
    pub fn from_voxels(density: Real, voxels: &Voxels) -> Self {
        let mut com = Point::origin();
        let mut num_not_empty = 0;
        let mut angular_inertia = na::zero();
        let block_ref_mprops = MassProperties::from_cuboid(density, voxels.voxel_size() / 2.0);

        for vox in voxels.voxels() {
            if !vox.state.is_empty() {
                com += vox.center.coords;
                num_not_empty += 1;
            }
        }

        com.coords /= num_not_empty as Real;

        for vox in voxels.voxels() {
            if !vox.state.is_empty() {
                angular_inertia +=
                    block_ref_mprops.construct_shifted_inertia_matrix(vox.center - com);
            }
        }

        let mass = block_ref_mprops.mass() * num_not_empty as Real;

        #[cfg(feature = "dim2")]
        return Self::new(com, mass, angular_inertia);
        #[cfg(feature = "dim3")]
        return Self::with_inertia_matrix(com, mass, angular_inertia);
    }
}
