use crate::mass_properties::MassProperties;
use crate::math::{Point, Real, DIM};

impl MassProperties {
    /// Computes the mass properties of a convex polyhedron.
    pub fn from_convex_polyhedron(
        density: Real,
        vertices: &[Point<Real>],
        indices: &[[u32; DIM]],
    ) -> MassProperties {
        Self::from_trimesh(density, vertices, indices)
    }
}
