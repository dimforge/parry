use crate::math::{Point, Real, Vector, DIM};
use crate::shape::Ball;
use crate::transformation::utils;
use na::{self, ComplexField, Point3, RealField};

impl Ball {
    /// Discretize the boundary of this ball as a triangle-mesh.
    pub fn to_trimesh(
        &self,
        ntheta_subdiv: u32,
        nphi_subdiv: u32,
    ) -> (Vec<Point3<Real>>, Vec<[u32; 3]>) {
        let diameter = self.radius * 2.0;
        let (vtx, idx) = unit_sphere(ntheta_subdiv, nphi_subdiv);
        (utils::scaled(vtx, Vector::repeat(diameter)), idx)
    }
}

fn unit_sphere(ntheta_subdiv: u32, nphi_subdiv: u32) -> (Vec<Point3<Real>>, Vec<[u32; 3]>) {
    let dtheta = Real::two_pi() / (ntheta_subdiv as Real);
    let dphi = Real::pi() / (nphi_subdiv as Real);

    let mut coords = Vec::new();
    let mut curr_phi: Real = -Real::frac_pi_2() + dphi;

    coords.push(Point::new(0.0, -1.0, 0.0));

    for _ in 1..nphi_subdiv {
        utils::push_circle(
            ComplexField::cos(curr_phi),
            ntheta_subdiv,
            dtheta,
            ComplexField::sin(curr_phi),
            &mut coords,
        );
        curr_phi += dphi;
    }

    coords.push(Point::new(0.0, 1.0, 0.0));

    let mut idx = Vec::new();

    utils::push_degenerate_top_ring_indices(1, 0, ntheta_subdiv, &mut idx);
    utils::reverse_clockwising(&mut idx);

    for i in 0..nphi_subdiv - 2 {
        utils::push_ring_indices(
            1 + i * ntheta_subdiv,
            1 + (i + 1) * ntheta_subdiv,
            ntheta_subdiv,
            &mut idx,
        );
    }

    utils::push_degenerate_top_ring_indices(
        coords.len() as u32 - 1 - ntheta_subdiv,
        coords.len() as u32 - 1,
        ntheta_subdiv,
        &mut idx,
    );

    (utils::scaled(coords, Vector::repeat(0.5)), idx)
}

/// Creates an hemisphere with a diameter of 1.
pub(crate) fn unit_hemisphere(
    ntheta_subdiv: u32,
    nphi_subdiv: u32,
) -> (Vec<Point<Real>>, Vec<[u32; DIM]>) {
    let two_pi = Real::two_pi();
    let pi_two = Real::frac_pi_2();
    let dtheta = two_pi / (ntheta_subdiv as Real);
    let dphi = pi_two / (nphi_subdiv as Real);

    let mut coords = Vec::new();
    let mut curr_phi: Real = 0.0;

    for _ in 0..nphi_subdiv {
        utils::push_circle(
            ComplexField::cos(curr_phi),
            ntheta_subdiv,
            dtheta,
            ComplexField::sin(curr_phi),
            &mut coords,
        );
        curr_phi += dphi;
    }

    coords.push(Point::new(0.0, 1.0, 0.0));

    let mut idx = Vec::new();

    for i in 0..nphi_subdiv - 1 {
        utils::push_ring_indices(
            i * ntheta_subdiv,
            (i + 1) * ntheta_subdiv,
            ntheta_subdiv,
            &mut idx,
        );
    }

    utils::push_degenerate_top_ring_indices(
        (nphi_subdiv - 1) * ntheta_subdiv,
        coords.len() as u32 - 1,
        ntheta_subdiv,
        &mut idx,
    );

    (utils::scaled(coords, Vector::repeat(0.5)), idx)
}
