use crate::math::{Isometry, Real, Vector, DIM};
use crate::shape::{Cuboid, SupportMap};

/// Computes the separation of two cuboids along `axis1`.
#[cfg(feature = "dim3")]
pub fn cuboid_cuboid_compute_separation_wrt_local_line(
    cuboid1: &Cuboid,
    cuboid2: &Cuboid,
    pos12: &Isometry<Real>,
    axis1: &Vector<Real>,
) -> (Real, Vector<Real>) {
    let signum = (1.0 as Real).copysign(pos12.translation.vector.dot(axis1));
    let axis1 = axis1 * signum;
    let axis2 = pos12.inverse_transform_vector(&-axis1);
    let local_pt1 = cuboid1.local_support_point(&axis1);
    let local_pt2 = cuboid2.local_support_point(&axis2);
    let pt2 = pos12 * local_pt2;
    let separation = (pt2 - local_pt1).dot(&axis1);
    (separation, axis1)
}

/// Finds the best separating edge between two cuboids.
///
/// All combinations of edges from both cuboids are taken into
/// account.
#[cfg(feature = "dim3")]
pub fn cuboid_cuboid_find_local_separating_edge_twoway(
    cuboid1: &Cuboid,
    cuboid2: &Cuboid,
    pos12: &Isometry<Real>,
) -> (Real, Vector<Real>) {
    use approx::AbsDiffEq;
    let mut best_separation = -Real::MAX;
    let mut best_dir = Vector::zeros();

    let x2 = pos12 * Vector::x();
    let y2 = pos12 * Vector::y();
    let z2 = pos12 * Vector::z();

    // We have 3 * 3 = 9 axes to test.
    let axes = [
        // Vector::{x, y ,z}().cross(y2)
        Vector::new(0.0, -x2.z, x2.y),
        Vector::new(x2.z, 0.0, -x2.x),
        Vector::new(-x2.y, x2.x, 0.0),
        // Vector::{x, y ,z}().cross(y2)
        Vector::new(0.0, -y2.z, y2.y),
        Vector::new(y2.z, 0.0, -y2.x),
        Vector::new(-y2.y, y2.x, 0.0),
        // Vector::{x, y ,z}().cross(y2)
        Vector::new(0.0, -z2.z, z2.y),
        Vector::new(z2.z, 0.0, -z2.x),
        Vector::new(-z2.y, z2.x, 0.0),
    ];

    for axis1 in &axes {
        let norm1 = axis1.norm();
        if norm1 > Real::default_epsilon() {
            let (separation, axis1) = cuboid_cuboid_compute_separation_wrt_local_line(
                cuboid1,
                cuboid2,
                pos12,
                &(axis1 / norm1),
            );

            if separation > best_separation {
                best_separation = separation;
                best_dir = axis1;
            }
        }
    }

    (best_separation, best_dir)
}

/// Finds the best separating normal between two cuboids.
///
/// Only the normals from `cuboid1` are tested.
pub fn cuboid_cuboid_find_local_separating_normal_oneway(
    cuboid1: &Cuboid,
    cuboid2: &Cuboid,
    pos12: &Isometry<Real>,
) -> (Real, Vector<Real>) {
    let mut best_separation = -Real::MAX;
    let mut best_dir = Vector::zeros();

    for i in 0..DIM {
        let sign = (1.0 as Real).copysign(pos12.translation.vector[i]);
        let axis1 = Vector::ith(i, sign);
        let axis2 = pos12.inverse_transform_vector(&-axis1);
        let local_pt2 = cuboid2.local_support_point(&axis2);
        let pt2 = pos12 * local_pt2;
        let separation = pt2[i] * sign - cuboid1.half_extents[i];

        if separation > best_separation {
            best_separation = separation;
            best_dir = axis1;
        }
    }

    (best_separation, best_dir)
}
