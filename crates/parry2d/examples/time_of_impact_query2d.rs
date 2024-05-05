extern crate nalgebra as na;

use na::{Isometry2, Vector2};
use parry2d::query;
use parry2d::query::ShapeCastOptions;
use parry2d::shape::{Ball, Cuboid};

fn main() {
    let cuboid = Cuboid::new(Vector2::new(1.0, 1.0));
    let ball = Ball::new(1.0);

    let cuboid_pos = Isometry2::identity();
    let ball_pos_intersecting = Isometry2::translation(1.0, 1.0);
    let ball_pos_will_touch = Isometry2::translation(2.0, 2.0);
    let ball_pos_wont_touch = Isometry2::translation(3.0, 3.0);

    let box_vel1 = Vector2::new(-1.0, 1.0);
    let box_vel2 = Vector2::new(1.0, 1.0);

    let ball_vel1 = Vector2::new(2.0, 2.0);
    let ball_vel2 = Vector2::new(-0.5, -0.5);

    let toi_intersecting = query::cast_shapes(
        &ball_pos_intersecting,
        &ball_vel1,
        &ball,
        &cuboid_pos,
        &box_vel1,
        &cuboid,
        ShapeCastOptions::default(),
    )
    .unwrap();
    let toi_will_touch = query::cast_shapes(
        &ball_pos_will_touch,
        &ball_vel2,
        &ball,
        &cuboid_pos,
        &box_vel2,
        &cuboid,
        ShapeCastOptions::default(),
    )
    .unwrap();
    let toi_wont_touch = query::cast_shapes(
        &ball_pos_wont_touch,
        &ball_vel1,
        &ball,
        &cuboid_pos,
        &box_vel1,
        &cuboid,
        ShapeCastOptions::default(),
    )
    .unwrap();

    assert_eq!(
        toi_intersecting.map(|time_of_impact| time_of_impact.time_of_impact),
        Some(0.0)
    );
    println!("Toi: {:?}", toi_will_touch);
    assert!(toi_will_touch.is_some() && toi_will_touch.unwrap().time_of_impact > 0.0);
    assert_eq!(
        toi_wont_touch.map(|time_of_impact| time_of_impact.time_of_impact),
        None
    );
}
