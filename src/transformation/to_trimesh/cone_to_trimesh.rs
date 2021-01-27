use crate::math::Real;
use crate::shape::Cone;
use crate::transformation::utils;
use na::{self, Point3, RealField, Vector3};

impl Cone {
    /// Discretize the boundary of this cone as a triangle-mesh.
    pub fn to_trimesh(&self, nsubdiv: u32) -> (Vec<Point3<Real>>, Vec<[u32; 3]>) {
        let diameter = self.radius * 2.0;
        let height = self.half_height * 2.0;
        let scale = Vector3::new(diameter, height, diameter);
        let (vtx, idx) = unit_cone(nsubdiv);
        (utils::scaled(vtx, scale), idx)
    }
}

/// Generates a cone with unit height and diameter.
fn unit_cone(nsubdiv: u32) -> (Vec<Point3<Real>>, Vec<[u32; 3]>) {
    let two_pi = Real::two_pi();
    let dtheta = two_pi / (nsubdiv as Real);
    let mut coords = Vec::new();
    let mut indices = Vec::new();

    utils::push_circle(
        na::convert(0.5),
        nsubdiv,
        dtheta,
        na::convert(-0.5),
        &mut coords,
    );

    coords.push(Point3::new(0.0, 0.5, 0.0));

    utils::push_degenerate_top_ring_indices(0, coords.len() as u32 - 1, nsubdiv, &mut indices);
    utils::push_filled_circle_indices(0, nsubdiv, &mut indices);

    (coords, indices)
}
