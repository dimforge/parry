use crate::math::Real;
use crate::shape::Cylinder;
use crate::transformation::utils;
use na::{self, Point3, RealField, Vector3};

impl Cylinder {
    /// Discretize the boundary of this cylinder as a triangle-mesh.
    pub fn to_outline(&self, nsubdiv: u32) -> (Vec<Point3<Real>>, Vec<[u32; 2]>) {
        let diameter = self.radius * 2.0;
        let height = self.half_height * 2.0;
        let scale = Vector3::new(diameter, height, diameter);
        let (vtx, idx) = unit_cylinder_outline(nsubdiv);
        (utils::scaled(vtx, scale), idx)
    }
}

/// Generates a cylinder with unit height and diameter.
fn unit_cylinder_outline(nsubdiv: u32) -> (Vec<Point3<Real>>, Vec<[u32; 2]>) {
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
    utils::push_circle(
        na::convert(0.5),
        nsubdiv,
        dtheta,
        na::convert(0.5),
        &mut coords,
    );

    utils::push_circle_outline_indices(&mut indices, 0..nsubdiv);
    utils::push_circle_outline_indices(&mut indices, nsubdiv..nsubdiv * 2);
    indices.push([0, nsubdiv]);
    indices.push([nsubdiv / 4, nsubdiv + nsubdiv / 4]);
    indices.push([nsubdiv / 2, nsubdiv + nsubdiv / 2]);
    indices.push([(nsubdiv * 3) / 4, nsubdiv + (nsubdiv * 3) / 4]);

    (coords, indices)
}
