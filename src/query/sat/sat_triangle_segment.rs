use crate::math::{Isometry, Real, Vector};
use crate::query::sat;
use crate::shape::{Segment, SupportMap, Triangle};
use na::Unit;

/// Finds the best separating normal a triangle and a segment.
///
/// Only the normals of `triangle1` are tested.
pub fn triangle_segment_find_local_separating_normal_oneway(
    triangle1: &Triangle,
    segment2: &Segment,
    pos12: &Isometry<Real>,
) -> (Real, Vector<Real>) {
    if let Some(dir) = triangle1.normal() {
        let p2a = segment2.support_point_toward(pos12, &-dir);
        let p2b = segment2.support_point_toward(pos12, &dir);
        let sep_a = (p2a - triangle1.a).dot(&dir);
        let sep_b = -(p2b - triangle1.a).dot(&dir);

        if sep_a >= sep_b {
            (sep_a, *dir)
        } else {
            (sep_b, -*dir)
        }
    } else {
        (-Real::MAX, Vector::zeros())
    }
}

/// Finds the best separating edge between a segment and a triangle.
///
/// All combinations of edges from the segment and the triangle are taken into
/// account.
pub fn segment_triangle_find_local_separating_edge_twoway(
    segment1: &Segment,
    triangle2: &Triangle,
    pos12: &Isometry<Real>,
) -> (Real, Vector<Real>) {
    let x2 = pos12 * (triangle2.b - triangle2.a);
    let y2 = pos12 * (triangle2.c - triangle2.b);
    let z2 = pos12 * (triangle2.a - triangle2.c);
    let dir1 = segment1.scaled_direction();

    let crosses1 = [dir1.cross(&x2), dir1.cross(&y2), dir1.cross(&z2)];
    let axes1 = [
        crosses1[0],
        crosses1[1],
        crosses1[2],
        -crosses1[0],
        -crosses1[1],
        -crosses1[2],
    ];
    let mut max_separation = -Real::MAX;
    let mut sep_dir = axes1[0];

    for axis1 in &axes1 {
        if let Some(axis1) = Unit::try_new(*axis1, 0.0) {
            let sep =
                sat::support_map_support_map_compute_separation(segment1, triangle2, pos12, &axis1);

            if sep > max_separation {
                max_separation = sep;
                sep_dir = *axis1;
            }
        }
    }

    (max_separation, sep_dir)
}
