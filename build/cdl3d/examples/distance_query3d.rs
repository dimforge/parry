#[macro_use]
extern crate approx; // for relative_eq!
extern crate nalgebra as na;

use cdl3d::query;
use cdl3d::shape::{Ball, Cuboid};
use na::{Isometry3, Vector3};

fn main() {
    let cuboid = Cuboid::new(Vector3::new(1.0, 1.0, 1.0));
    let ball = Ball::new(1.0);

    let cuboid_pos = Isometry3::identity();
    let ball_pos_intersecting = Isometry3::translation(0.0, 1.0, 0.0);
    let ball_pos_disjoint = Isometry3::translation(0.0, 3.0, 0.0);

    let dist_intersecting =
        query::distance(&ball_pos_intersecting, &ball, &cuboid_pos, &cuboid).unwrap();
    let dist_disjoint = query::distance(&ball_pos_disjoint, &ball, &cuboid_pos, &cuboid).unwrap();

    assert_eq!(dist_intersecting, 0.0);
    assert!(relative_eq!(dist_disjoint, 1.0, epsilon = 1.0e-7));
}
