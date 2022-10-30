use crate::bounding_volume::Aabb;
use crate::math::{Point, Real, Vector};
use crate::shape::Cuboid;
use crate::transformation::utils;

impl Aabb {
    /// Outlines this Aabb’s shape using polylines.
    pub fn to_outline(&self) -> (Vec<Point<Real>>, Vec<[u32; 2]>) {
        let center = self.center();
        let half_extents = self.half_extents();
        let mut cube_mesh = Cuboid::new(half_extents).to_outline();
        cube_mesh.0.iter_mut().for_each(|p| *p += center.coords);
        cube_mesh
    }
}

impl Cuboid {
    /// Outlines this cuboid’s shape using polylines.
    pub fn to_outline(&self) -> (Vec<Point<Real>>, Vec<[u32; 2]>) {
        let (vtx, idx) = unit_cuboid_outline();
        (utils::scaled(vtx, self.half_extents * 2.0), idx)
    }
}

/**
 * Generates a cuboid shape with a split index buffer.
 *
 * The cuboid is centered at the origin, and has its half extents set to 0.5.
 */
fn unit_cuboid_outline() -> (Vec<Point<Real>>, Vec<[u32; 2]>) {
    let aabb = Aabb::from_half_extents(Point::origin(), Vector::repeat(0.5));
    (
        aabb.vertices().to_vec(),
        vec![
            [0, 1],
            [1, 2],
            [2, 3],
            [3, 0],
            [4, 5],
            [5, 6],
            [6, 7],
            [7, 4],
            [0, 4],
            [1, 5],
            [2, 6],
            [3, 7],
        ],
    )
}
