use parry3d::math::{Pose, Real, Vector};
use parry3d::query::{self, ShapeCastOptions};
use parry3d::shape::{Ball, Cuboid};

#[test]
fn ball_cuboid_toi() {
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
    .unwrap()
    .map(|hit| hit.time_of_impact);
    let toi_will_touch = query::cast_shapes(
        &ball_pos_will_touch,
        ball_vel2,
        &ball,
        &cuboid_pos,
        cuboid_vel2,
        &cuboid,
        ShapeCastOptions::default(),
    )
    .unwrap()
    .map(|hit| hit.time_of_impact);
    let toi_wont_touch = query::cast_shapes(
        &ball_pos_wont_touch,
        ball_vel1,
        &ball,
        &cuboid_pos,
        cuboid_vel1,
        &cuboid,
        ShapeCastOptions::default(),
    )
    .unwrap()
    .map(|hit| hit.time_of_impact);

    assert_eq!(toi_intersecting, Some(0.0));
    assert!(relative_eq!(
        toi_will_touch.unwrap(),
        ((3.0 as Real).sqrt() - 1.0) / (ball_vel2 - cuboid_vel2).length()
    ));
    assert_eq!(toi_wont_touch, None);
}
