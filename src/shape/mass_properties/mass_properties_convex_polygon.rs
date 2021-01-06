#![allow(dead_code)] // TODO: remove this

use crate::math::{Point, Real};
use crate::shape::{MassProperties, Triangle};

impl MassProperties {
    pub(crate) fn from_convex_polygon(density: Real, vertices: &[Point<Real>]) -> MassProperties {
        let (area, com) = convex_polygon_area_and_center_of_mass(vertices);

        if area == 0.0 {
            return MassProperties::new(com, 0.0, 0.0);
        }

        let mut itot = 0.0;
        let factor = 1.0 / 6.0;

        let mut iterpeek = vertices.iter().peekable();
        let firstelement = *iterpeek.peek().unwrap(); // store first element to close the cycle in the end with unwrap_or
        while let Some(elem) = iterpeek.next() {
            let triangle = Triangle::new(com, *elem, **iterpeek.peek().unwrap_or(&firstelement));
            let area = triangle.area();

            // algorithm adapted from Box2D
            let e1 = triangle.b - com;
            let e2 = triangle.c - com;

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

fn convex_polygon_area_and_center_of_mass(convex_polygon: &[Point<Real>]) -> (Real, Point<Real>) {
    let geometric_center = convex_polygon
        .iter()
        .fold(Point::origin(), |e1, e2| e1 + e2.coords)
        / convex_polygon.len() as Real;
    let mut res = Point::origin();
    let mut areasum = 0.0;

    let mut iterpeek = convex_polygon.iter().peekable();
    let firstelement = *iterpeek.peek().unwrap(); // Stores first element to close the cycle in the end with unwrap_or.
    while let Some(elem) = iterpeek.next() {
        let (a, b, c) = (
            elem,
            iterpeek.peek().unwrap_or(&firstelement),
            &geometric_center,
        );
        let area = Triangle::new(*a, **b, *c).area();
        let center = (a.coords + b.coords + c.coords) / 3.0;

        res += center * area;
        areasum += area;
    }

    if areasum == 0.0 {
        (areasum, geometric_center)
    } else {
        (areasum, res / areasum)
    }
}
