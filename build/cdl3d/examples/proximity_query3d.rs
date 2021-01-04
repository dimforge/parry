extern crate nalgebra as na;

use cdl3d::query;
use cdl3d::shape::{Ball, Cuboid};
use na::{Isometry3, Vector3};

fn main() {
    let cuboid = Cuboid::new(Vector3::new(1.0, 1.0, 1.0));
    let ball = Ball::new(1.0);
    let cuboid_pos = Isometry3::identity();
    let ball_pos_intersecting = Isometry3::translation(1.0, 1.0, 1.0);
    let ball_pos_disjoint = Isometry3::translation(3.0, 3.0, 3.0);

    let intersecting =
        query::intersection_test(&ball_pos_intersecting, &ball, &cuboid_pos, &cuboid).unwrap();
    let not_intersecting =
        !query::intersection_test(&ball_pos_disjoint, &ball, &cuboid_pos, &cuboid).unwrap();

    assert!(intersecting);
    assert!(not_intersecting);
}
