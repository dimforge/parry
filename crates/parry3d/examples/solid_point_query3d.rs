extern crate nalgebra as na;

use parry3d::math::{Isometry, Point, PointOps, Vector};
use parry3d::query::PointQuery;
use parry3d::shape::Cuboid;

fn main() {
    let cuboid = Cuboid::new(Vector::new(1.0, 2.0, 2.0));
    let pt_inside = Point::origin();
    let pt_outside = Point::new(2.0, 2.0, 2.0);

    // Solid projection.
    assert_eq!(
        cuboid.distance_to_point(&Isometry::identity(), &pt_inside, true),
        0.0
    );

    // Non-solid projection.
    assert_eq!(
        cuboid.distance_to_point(&Isometry::identity(), &pt_inside, false),
        -1.0
    );

    // The other point is outside of the cuboid so the `solid` flag has no effect.
    assert_eq!(
        cuboid.distance_to_point(&Isometry::identity(), &pt_outside, false),
        1.0
    );
    assert_eq!(
        cuboid.distance_to_point(&Isometry::identity(), &pt_outside, true),
        1.0
    );
}
