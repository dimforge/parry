#![allow(dead_code)] // TODO: remove this

use crate::math::{Point, Real};
use crate::shape::{MassProperties, Triangle};
use na::Point3;

impl MassProperties {
    pub(crate) fn from_trimesh(
        density: Real,
        vertices: &[Point<Real>],
        indices: &[Point3<u32>],
    ) -> MassProperties {
        let (area, com) = trimesh_area_and_center_of_mass(vertices, indices);

        if area == 0.0 {
            return MassProperties::new(com, 0.0, 0.0);
        }

        let mut itot = 0.0;
        let factor = 1.0 / 6.0;

        for idx in indices {
            let triangle = Triangle::new(
                vertices[idx.x as usize],
                vertices[idx.y as usize],
                vertices[idx.z as usize],
            );
            let area = triangle.area();

            // TODO: is the parallel axis theorem correctly applied here?
            // algorithm adapted from the convex polygon code.
            let e1 = triangle.b - triangle.a;
            let e2 = triangle.c - triangle.a;

            let ex1 = e1[0];
            let ey1 = e1[1];
            let ex2 = e2[0];
            let ey2 = e2[1];

            let intx2 = ex1 * ex1 + ex2 * ex1 + ex2 * ex2;
            let inty2 = ey1 * ey1 + ey2 * ey1 + ey2 * ey2;

            let ipart = factor * (intx2 + inty2);

            itot += ipart * area;
        }

        Self::new(com, area * density, itot * density)
    }
}

fn trimesh_area_and_center_of_mass(
    vertices: &[Point<Real>],
    indices: &[Point3<u32>],
) -> (Real, Point<Real>) {
    let mut res = Point::origin();
    let mut areasum = 0.0;

    for idx in indices {
        let triangle = Triangle::new(
            vertices[idx.x as usize],
            vertices[idx.y as usize],
            vertices[idx.z as usize],
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
