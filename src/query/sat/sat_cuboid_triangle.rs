use crate::math::{Isometry, Real, Vector};
use crate::query::sat;
use crate::query::sat::support_map_support_map_compute_separation;
use crate::shape::{Cuboid, SupportMap, Triangle};

#[cfg(feature = "dim3")]
pub fn cuboid_triangle_find_local_separating_edge_twoway(
    cube1: &Cuboid,
    triangle2: &Triangle,
    pos12: &Isometry<Real>,
) -> (Real, Vector<Real>) {
    let x2 = pos12 * (triangle2.b - triangle2.a);
    let y2 = pos12 * (triangle2.c - triangle2.b);
    let z2 = pos12 * (triangle2.a - triangle2.c);

    // We have 3 * 3 = 3 axes to test.
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

    sat::cuboid_support_map_find_local_separating_edge_twoway(cube1, triangle2, &axes, pos12)
}

#[cfg(feature = "dim2")]
pub fn triangle_support_map_find_local_separating_normal_oneway(
    triangle1: &Triangle,
    shape2: &impl SupportMap,
    pos12: &Isometry<Real>,
) -> (Real, Vector<Real>) {
    let mut best_sep = -Real::MAX;
    let mut best_normal = Vector::zeros();

    for edge in &triangle1.edges() {
        if let Some(normal) = edge.normal() {
            let sep = support_map_support_map_compute_separation(triangle1, shape2, pos12, &normal);

            if sep > best_sep {
                best_sep = sep;
                best_normal = *normal;
            }
        }
    }

    (best_sep, best_normal)
}

#[cfg(feature = "dim2")]
pub fn triangle_cuboid_find_local_separating_normal_oneway(
    triangle1: &Triangle,
    shape2: &Cuboid,
    pos12: &Isometry<Real>,
) -> (Real, Vector<Real>) {
    triangle_support_map_find_local_separating_normal_oneway(triangle1, shape2, pos12)
}

#[cfg(feature = "dim3")]
pub fn triangle_cuboid_find_local_separating_normal_oneway(
    triangle1: &Triangle,
    shape2: &Cuboid,
    pos12: &Isometry<Real>,
) -> (Real, Vector<Real>) {
    sat::point_cuboid_find_local_separating_normal_oneway(
        triangle1.a,
        triangle1.normal(),
        shape2,
        pos12,
    )
}
