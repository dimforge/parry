use crate::bounding_volume::AABB;
use crate::math::{Isometry, Real, Translation};
use crate::query::sat;
use crate::shape::{Cuboid, Triangle};

// TODO: put this somewhere in the `query` module?
pub fn aabb_intersects_triangle(aabb1: &AABB, triangle2: &Triangle) -> bool {
    let cube1 = Cuboid::new(aabb1.half_extents());
    let pos12 = Isometry::from(Translation::from(-aabb1.center().coords));
    let pos21 = Isometry::from(Translation::from(aabb1.center().coords));

    /*
     *
     * Point-Face cases.
     *
     */
    let sep1 =
        sat::cuboid_support_map_find_local_separating_normal_oneway(&cube1, triangle2, &pos12).0;
    if sep1 > 0.0 {
        return false;
    }

    let sep2 =
        sat::triangle_cuboid_find_local_separating_normal_oneway(triangle2, &cube1, &pos21).0;
    if sep2 > 0.0 {
        return false;
    }

    /*
     *
     * Edge-Edge cases.
     *
     */
    #[cfg(feature = "dim2")]
    let sep3 = -Real::MAX; // This case does not exist in 2D.
    #[cfg(feature = "dim3")]
    let sep3 = sat::cuboid_triangle_find_local_separating_edge_twoway(&cube1, triangle2, &pos12).0;

    sep3 <= 0.0
}
