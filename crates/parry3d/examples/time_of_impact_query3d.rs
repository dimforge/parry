use parry3d::math::{Pose, Vector};
use parry3d::query::{self, ShapeCastOptions};
use parry3d::shape::{Ball, Cuboid};

fn main() {
    let cuboid = Cuboid::new(Vector::new(1.0, 1.0, 1.0));
    let ball = Ball::new(1.0);

    let cuboid_pos = Pose::identity();
    let ball_pos_intersecting = Pose::translation(1.0, 1.0, 1.0);
    let ball_pos_will_touch = Pose::translation(2.0, 2.0, 2.0);
    let ball_pos_wont_touch = Pose::translation(3.0, 3.0, 3.0);

    let cuboid_vel1 = Vector::new(-1.0, 1.0, 1.0);
    let cuboid_vel2 = Vector::new(1.0, 1.0, 1.0);

    let ball_vel1 = Vector::new(2.0, 2.0, 2.0);
    let ball_vel2 = Vector::new(-0.5, -0.5, -0.5);

    let toi_intersecting = query::cast_shapes(
        &ball_pos_intersecting,
        ball_vel1,
        &ball,
        &cuboid_pos,
        cuboid_vel1,
        &cuboid,
        ShapeCastOptions::default(),
    )
    .unwrap();
    let toi_will_touch = query::cast_shapes(
        &ball_pos_will_touch,
        ball_vel2,
        &ball,
        &cuboid_pos,
        cuboid_vel2,
        &cuboid,
        ShapeCastOptions::default(),
    )
    .unwrap();
    let toi_wont_touch = query::cast_shapes(
        &ball_pos_wont_touch,
        ball_vel1,
        &ball,
        &cuboid_pos,
        cuboid_vel1,
        &cuboid,
        ShapeCastOptions::default(),
    )
    .unwrap();

    assert_eq!(toi_intersecting.map(|hit| hit.time_of_impact), Some(0.0));
    assert!(toi_will_touch.is_some() && toi_will_touch.unwrap().time_of_impact > 0.0);
    assert_eq!(toi_wont_touch.map(|hit| hit.time_of_impact), None);
}
