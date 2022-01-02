extern crate nalgebra as na;

use na::{Isometry2, Vector2};
use parry2d::query;
use parry2d::shape::{Ball, Cuboid};

fn main() {
    let cuboid = Cuboid::new(Vector2::new(1.0, 1.0));
    let ball = Ball::new(1.0);

    let cuboid_pos = Isometry2::identity();
    let ball_pos_intersecting = Isometry2::translation(1.0, 1.0);
    let ball_pos_disjoint = Isometry2::translation(3.0, 3.0);

    assert!(query::intersection_test(&ball_pos_intersecting, &ball, &cuboid_pos, &cuboid).unwrap());
    assert!(!query::intersection_test(&ball_pos_disjoint, &ball, &cuboid_pos, &cuboid).unwrap());
}
