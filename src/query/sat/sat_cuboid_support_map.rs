use crate::math::{Isometry, Real, Vector, DIM};
use crate::shape::{Cuboid, SupportMap};

use na::Unit;

/// Computes the separation between a cuboid an a convex shape implementing the `SupportMap` trait,
/// along the given axis.
// TODO: this is a very slow approach. We should only do special-cases instead.
#[cfg(feature = "dim3")]
pub fn cuboid_support_map_compute_separation_wrt_local_line(
    cube1: &Cuboid,
    shape2: &impl SupportMap,
    pos12: &Isometry<Real>,
    axis1: &Unit<Vector<Real>>,
) -> (Real, Unit<Vector<Real>>) {
    let axis1_2 = pos12.inverse_transform_unit_vector(&axis1);
    let separation1 = {
        let axis2 = -axis1_2;
        let local_pt1 = cube1.local_support_point_toward(&axis1);
        let local_pt2 = shape2.local_support_point_toward(&axis2);
        let pt2 = pos12 * local_pt2;
        (pt2 - local_pt1).dot(&axis1)
    };

    let separation2 = {
        let axis2 = axis1_2;
        let local_pt1 = cube1.local_support_point_toward(&-*axis1);
        let local_pt2 = shape2.local_support_point_toward(&axis2);
        let pt2 = pos12 * local_pt2;
        (pt2 - local_pt1).dot(&-*axis1)
    };

    if separation1 > separation2 {
        (separation1, *axis1)
    } else {
        (separation2, -*axis1)
    }
}

/// Finds the best separating edge between a cuboid and a convex shape implementing the `Supportmap` trait.
///
/// Only the axes given by `axes` are tested.
#[cfg(feature = "dim3")]
pub fn cuboid_support_map_find_local_separating_edge_twoway(
    cube1: &Cuboid,
    shape2: &impl SupportMap,
    axes: &[Vector<Real>],
    pos12: &Isometry<Real>,
) -> (Real, Vector<Real>) {
    use approx::AbsDiffEq;
    let mut best_separation = -Real::MAX;
    let mut best_dir = Vector::zeros();

    for axis1 in axes {
        if let Some(axis1) = Unit::try_new(*axis1, Real::default_epsilon()) {
            let (separation, axis1) =
                cuboid_support_map_compute_separation_wrt_local_line(cube1, shape2, pos12, &axis1);

            if separation > best_separation {
                best_separation = separation;
                best_dir = *axis1;
            }
        }
    }

    (best_separation, best_dir)
}

/// Finds the best separating normal between a cuboid and a convex shape implementing the `SupportMap` trait.
///
/// Only the normals of `cube1` are tested.
pub fn cuboid_support_map_find_local_separating_normal_oneway<S: SupportMap>(
    cube1: &Cuboid,
    shape2: &S,
    pos12: &Isometry<Real>,
) -> (Real, Vector<Real>) {
    let mut best_separation = -Real::MAX;
    let mut best_dir = Vector::zeros();

    for i in 0..DIM {
        for sign in &[-1.0, 1.0] {
            let axis1 = Vector::ith(i, *sign);
            let pt2 = shape2.support_point_toward(&pos12, &Unit::new_unchecked(-axis1));
            let separation = pt2[i] * *sign - cube1.half_extents[i];

            if separation > best_separation {
                best_separation = separation;
                best_dir = axis1;
            }
        }
    }

    (best_separation, best_dir)
}
