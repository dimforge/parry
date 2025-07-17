extern crate nalgebra as na;

use na::{Isometry2, Point2, Vector2};
use parry2d::query::PointQuery;
use parry2d::shape::Cuboid;

fn main() {
    let cuboid = Cuboid::new(Vector2::new(1.0, 2.0));
    let pt_inside = Point2::origin();
    let pt_outside = Point2::new(2.0, 2.0);

    let options = &();
    // Solid projection.
    assert_eq!(
        cuboid.distance_to_point(&Isometry2::identity(), &pt_inside, true, options),
        0.0
    );

    // Non-solid projection.
    assert_eq!(
        cuboid.distance_to_point(&Isometry2::identity(), &pt_inside, false, options),
        -1.0
    );

    // The other point is outside of the cuboid so the `solid` flag has no effect.
    assert_eq!(
        cuboid.distance_to_point(&Isometry2::identity(), &pt_outside, false, options),
        1.0
    );
    assert_eq!(
        cuboid.distance_to_point(&Isometry2::identity(), &pt_outside, true, options),
        1.0
    );
}
