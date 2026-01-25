use parry2d::math::{Pose, Vector};
use parry2d::query;
use parry2d::shape::{Ball, Cuboid};

fn main() {
    let cuboid = Cuboid::new(Vector::new(1.0, 1.0));
    let ball = Ball::new(1.0);

    let cuboid_pos = Pose::identity();
    let ball_pos_intersecting = Pose::translation(1.0, 1.0);
    let ball_pos_disjoint = Pose::translation(3.0, 3.0);

    assert!(query::intersection_test(&ball_pos_intersecting, &ball, &cuboid_pos, &cuboid).unwrap());
    assert!(!query::intersection_test(&ball_pos_disjoint, &ball, &cuboid_pos, &cuboid).unwrap());
}
