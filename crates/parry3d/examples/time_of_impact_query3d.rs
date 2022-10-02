extern crate nalgebra as na;

use na::{Isometry3, Vector3};
use parry3d::math::Real;
use parry3d::query;
use parry3d::shape::{Ball, Cuboid};

fn main() {
    let cuboid = Cuboid::new(Vector3::new(1.0, 1.0, 1.0));
    let ball = Ball::new(1.0);

    let cuboid_pos = Isometry3::identity();
    let ball_pos_intersecting = Isometry3::translation(1.0, 1.0, 1.0);
    let ball_pos_will_touch = Isometry3::translation(2.0, 2.0, 2.0);
    let ball_pos_wont_touch = Isometry3::translation(3.0, 3.0, 3.0);

    let cuboid_vel1 = Vector3::new(-1.0, 1.0, 1.0);
    let cuboid_vel2 = Vector3::new(1.0, 1.0, 1.0);

    let ball_vel1 = Vector3::new(2.0, 2.0, 2.0);
    let ball_vel2 = Vector3::new(-0.5, -0.5, -0.5);

    let toi_intersecting = query::time_of_impact(
        &ball_pos_intersecting,
        &ball_vel1,
        &ball,
        &cuboid_pos,
        &cuboid_vel1,
        &cuboid,
        Real::MAX,
        true,
    )
    .unwrap();
    let toi_will_touch = query::time_of_impact(
        &ball_pos_will_touch,
        &ball_vel2,
        &ball,
        &cuboid_pos,
        &cuboid_vel2,
        &cuboid,
        Real::MAX,
        true,
    )
    .unwrap();
    let toi_wont_touch = query::time_of_impact(
        &ball_pos_wont_touch,
        &ball_vel1,
        &ball,
        &cuboid_pos,
        &cuboid_vel1,
        &cuboid,
        Real::MAX,
        true,
    )
    .unwrap();

    assert_eq!(toi_intersecting.map(|toi| toi.toi), Some(0.0));
    assert!(toi_will_touch.is_some() && toi_will_touch.unwrap().toi > 0.0);
    assert_eq!(toi_wont_touch.map(|toi| toi.toi), None);
}
