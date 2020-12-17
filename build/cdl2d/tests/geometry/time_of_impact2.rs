use cdl2d::math::Real;
use cdl2d::query;
use cdl2d::shape::{Ball, Cuboid};
use na::{self, Isometry2, Vector2};

#[test]
fn ball_cuboid_toi() {
    let cuboid = Cuboid::new(Vector2::new(1.0, 1.0));
    let ball = Ball::new(1.0);

    let cuboid_pos = Isometry2::identity();
    let ball_pos_intersecting = Isometry2::new(Vector2::new(1.0, 1.0), 0.0);
    let ball_pos_will_touch = Isometry2::new(Vector2::new(2.0, 2.0), 0.0);
    let ball_pos_wont_touch = Isometry2::new(Vector2::new(3.0, 3.0), 0.0);

    let cuboid_vel1 = Vector2::new(-1.0, 1.0);
    let cuboid_vel2 = Vector2::new(1.0, 1.0);

    let ball_vel1 = Vector2::new(2.0, 2.0);
    let ball_vel2 = Vector2::new(-0.5, -0.5);

    let toi_intersecting = query::time_of_impact(
        &(ball_pos_intersecting.inverse() * cuboid_pos),
        &(cuboid_vel1 - ball_vel1),
        &ball,
        &cuboid,
        Real::MAX,
        0.0,
    )
    .unwrap();
    let toi_will_touch = query::time_of_impact(
        &(ball_pos_will_touch.inverse() * cuboid_pos),
        &(cuboid_vel2 - ball_vel2),
        &ball,
        &cuboid,
        Real::MAX,
        0.0,
    )
    .unwrap();
    let toi_wont_touch = query::time_of_impact(
        &(ball_pos_wont_touch.inverse() * cuboid_pos),
        &(cuboid_vel1 - ball_vel2),
        &ball,
        &cuboid,
        Real::MAX,
        0.0,
    )
    .unwrap();

    assert_eq!(toi_intersecting.map(|toi| toi.toi), Some(0.0));
    assert!(relative_eq!(
        toi_will_touch.unwrap().toi,
        ((2.0 as Real).sqrt() - 1.0) / (ball_vel2 - cuboid_vel2).norm()
    ));
    assert_eq!(toi_wont_touch.map(|toi| toi.toi), None);
}

#[test]
fn cuboid_cuboid_toi_issue_214() {
    let shape1 = Cuboid::new(Vector2::new(1.0, 1.0));
    let shape2 = Cuboid::new(Vector2::new(1.0, 1.5));

    let pos1 = Isometry2::new(Vector2::new(0.0, 0.0), 0.0);
    let pos2 = Isometry2::new(Vector2::new(10.0, 0.0), 0.0);
    let pos12 = pos1.inverse() * pos2;

    let vel1 = Vector2::new(1.0, 0.0);
    let vel2 = Vector2::new(0.0, 0.0);
    let vel12 = vel2 - vel1;

    let toi = query::time_of_impact(&pos12, &vel12, &shape1, &shape2, Real::MAX, 0.0).unwrap();
    assert!(toi.is_some());
}
