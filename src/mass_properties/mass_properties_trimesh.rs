use crate::mass_properties::MassProperties;
use crate::math::{Point, Real};
use crate::shape::Triangle;

impl MassProperties {
    /// Computes the mass properties of a triangle-mesh.
    pub fn from_trimesh(
        density: Real,
        vertices: &[Point<Real>],
        indices: &[[u32; 3]],
    ) -> MassProperties {
        let (area, com) = trimesh_area_and_center_of_mass(vertices, indices);

        if area == 0.0 {
            return MassProperties::new(com, 0.0, 0.0);
        }

        let mut itot = 0.0;

        for idx in indices {
            let triangle = Triangle::new(
                vertices[idx[0] as usize],
                vertices[idx[1] as usize],
                vertices[idx[2] as usize],
            );

            // TODO: is the parallel axis theorem correctly applied here?
            let area = triangle.area();
            let ipart = triangle.unit_angular_inertia();
            itot += ipart * area;
        }

        Self::new(com, area * density, itot * density)
    }
}

/// Computes the area and center-of-mass of a triangle-mesh.
pub fn trimesh_area_and_center_of_mass(
    vertices: &[Point<Real>],
    indices: &[[u32; 3]],
) -> (Real, Point<Real>) {
    let mut res = Point::origin();
    let mut areasum = 0.0;

    for idx in indices {
        let triangle = Triangle::new(
            vertices[idx[0] as usize],
            vertices[idx[1] as usize],
            vertices[idx[2] as usize],
        );
        let area = triangle.area();
        let center = triangle.center();

        res += center.coords * area;
        areasum += area;
    }

    if areasum == 0.0 {
        (areasum, res)
    } else {
        (areasum, res / areasum)
    }
}
