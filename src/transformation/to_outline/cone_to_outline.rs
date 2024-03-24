use crate::math::Real;
use crate::shape::Cone;
use crate::transformation::utils;
use na::{self, Point3, Vector3};

impl Cone {
    /// Outlines this cone’s shape using polylines.
    pub fn to_outline(&self, nsubdiv: u32) -> (Vec<Point3<Real>>, Vec<[u32; 2]>) {
        let diameter = self.radius * 2.0;
        let height = self.half_height * 2.0;
        let scale = Vector3::new(diameter, height, diameter);
        let (vtx, idx) = unit_cone_outline(nsubdiv);
        (utils::scaled(vtx, scale), idx)
    }
}

/// Generates a cone with unit height and diameter.
fn unit_cone_outline(nsubdiv: u32) -> (Vec<Point3<Real>>, Vec<[u32; 2]>) {
    let mut out_vtx = vec![Point3::new(-0.5, -0.5, 0.0), Point3::new(0.0, 0.5, 0.0)];
    let mut out_ptx = vec![];
    #[allow(clippy::single_range_in_vec_init)] // The single range is on purpose.
    utils::apply_revolution(false, true, &[0..1], nsubdiv, &mut out_vtx, &mut out_ptx);
    (out_vtx, out_ptx)
}
