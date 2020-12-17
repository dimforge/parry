use crate::math::{Isometry, Real, Vector};
use crate::query::sat;
use crate::shape::{Cuboid, Segment};

#[cfg(feature = "dim3")]
pub fn cuboid_segment_find_local_separating_edge_twoway(
    cube1: &Cuboid,
    segment2: &Segment,
    pos12: &Isometry<Real>,
) -> (Real, Vector<Real>) {
    let x2 = pos12 * (segment2.b - segment2.a);

    let axes = [
        // Vector::{x, y ,z}().cross(y2)
        Vector::new(0.0, -x2.z, x2.y),
        Vector::new(x2.z, 0.0, -x2.x),
        Vector::new(-x2.y, x2.x, 0.0),
    ];

    sat::cuboid_support_map_find_local_separating_edge_twoway(cube1, segment2, &axes, pos12)
}

#[cfg(feature = "dim2")]
pub fn segment_cuboid_find_local_separating_normal_oneway(
    segment1: &Segment,
    shape2: &Cuboid,
    pos12: &Isometry<Real>,
) -> (Real, Vector<Real>) {
    sat::point_cuboid_find_local_separating_normal_oneway(
        segment1.a,
        segment1.normal(),
        shape2,
        pos12,
    )
}
