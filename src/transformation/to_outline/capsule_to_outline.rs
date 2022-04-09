use crate::math::Real;
use crate::shape::{Capsule, Cylinder};
use crate::transformation::utils;
use na::{self, Point3};

impl Capsule {
    /// Discretize the boundary of this capsule as a triangle-mesh.
    pub fn to_outline(&self, nsubdiv: u32) -> (Vec<Point3<Real>>, Vec<[u32; 2]>) {
        let (vtx, idx) = canonical_capsule_outline(self.radius, self.half_height(), nsubdiv);
        (utils::transformed(vtx, self.canonical_transform()), idx)
    }
}

/// Generates a capsule.
pub(crate) fn canonical_capsule_outline(
    caps_radius: Real,
    cylinder_half_height: Real,
    nsubdiv: u32,
) -> (Vec<Point3<Real>>, Vec<[u32; 2]>) {
    let (mut vtx, mut idx) = Cylinder::new(cylinder_half_height, caps_radius).to_outline(nsubdiv);

    // Generate the hemispheres.
    super::ball_to_outline::push_unit_hemisphere_outline(nsubdiv / 2, &mut vtx, &mut idx);
    super::ball_to_outline::push_unit_hemisphere_outline(nsubdiv / 2, &mut vtx, &mut idx);

    let ncap_pts = (nsubdiv / 2 + 1) * 2;
    vtx[(nsubdiv * 2) as usize..(nsubdiv * 2 + ncap_pts) as usize]
        .iter_mut()
        .for_each(|pt| {
            *pt *= caps_radius * 2.0;
            pt.y += cylinder_half_height
        });

    vtx[(nsubdiv * 2 + ncap_pts) as usize..]
        .iter_mut()
        .for_each(|pt| {
            *pt *= caps_radius * 2.0;
            pt.y = -pt.y - cylinder_half_height
        });

    (vtx, idx)
}
