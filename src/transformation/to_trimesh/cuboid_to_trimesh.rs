use crate::math::Real;
use crate::shape::Cuboid;
use crate::transformation::utils;
use na::{self, Point3};

impl Cuboid {
    pub fn to_trimesh(&self) -> (Vec<Point3<Real>>, Vec<[u32; 3]>) {
        let (vtx, idx) = unit_cuboid();
        (utils::scaled(vtx, self.half_extents * 2.0), idx)
    }
}

/**
 * Generates a cuboid shape with a split index buffer.
 *
 * The cuboid is centered at the origin, and has its half extents set to 0.5.
 */
fn unit_cuboid() -> (Vec<Point3<Real>>, Vec<[u32; 3]>) {
    let mut coords = Vec::with_capacity(8);
    let mut faces = Vec::with_capacity(12);

    coords.push(Point3::new(-0.5, -0.5, 0.5));
    coords.push(Point3::new(-0.5, -0.5, -0.5));
    coords.push(Point3::new(0.5, -0.5, -0.5));
    coords.push(Point3::new(0.5, -0.5, 0.5));
    coords.push(Point3::new(-0.5, 0.5, 0.5));
    coords.push(Point3::new(-0.5, 0.5, -0.5));
    coords.push(Point3::new(0.5, 0.5, -0.5));
    coords.push(Point3::new(0.5, 0.5, 0.5));

    faces.push([4, 5, 0]);
    faces.push([5, 1, 0]);

    faces.push([5, 6, 1]);
    faces.push([6, 2, 1]);

    faces.push([6, 7, 3]);
    faces.push([2, 6, 3]);

    faces.push([7, 4, 0]);
    faces.push([3, 7, 0]);

    faces.push([0, 1, 2]);
    faces.push([3, 0, 2]);

    faces.push([7, 6, 5]);
    faces.push([4, 7, 5]);

    (coords, faces)
}
