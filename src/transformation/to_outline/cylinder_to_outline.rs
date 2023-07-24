use crate::math::{Real, real};
use crate::shape::Cylinder;
use crate::transformation::utils;
use na::{self, Point3, Vector3};

impl Cylinder {
    /// Outlines this cylinder’s shape using polylines.
    pub fn to_outline(&self, nsubdiv: u32) -> (Vec<Point3<Real>>, Vec<[u32; 2]>) {
        let diameter = self.radius * real!(2.0);
        let height = self.half_height * real!(2.0);
        let scale = Vector3::new(diameter, height, diameter);
        let (vtx, idx) = unit_cylinder_outline(nsubdiv);
        (utils::scaled(vtx, scale), idx)
    }
}

/// Generates a cylinder with unit height and diameter.
fn unit_cylinder_outline(nsubdiv: u32) -> (Vec<Point3<Real>>, Vec<[u32; 2]>) {
    let mut out_vtx = vec![Point3::new(real!(-0.5), real!(-0.5), real!(0.0)), Point3::new(real!(-0.5), real!(0.5), real!(0.0))];
    let mut out_idx = vec![];
    utils::apply_revolution(false, false, &[0..2], nsubdiv, &mut out_vtx, &mut out_idx);
    (out_vtx, out_idx)
}
