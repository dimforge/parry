use na::{self, Isometry2, Point2, Vector2};
use parry2d::math::Real;
use parry2d::query;
use parry2d::query::details::ShapeCastOptions;
use parry2d::shape::{Ball, Cuboid, Polyline, Segment};

#[test]
fn ball_ball_intersecting_toi() {
    let ball1 = Ball::new(1.0);
    let ball2 = Ball::new(2.0);

    let ball1_pos_intersecting = Isometry2::new(Vector2::new(1.0, 1.0), 0.0);
    let ball2_pos = Isometry2::identity();

    let ball1_vel_separating = Vector2::new(1.0, 2.0);
    let ball1_vel_penetrating = Vector2::new(1.0, -2.0);
    let ball2_vel = Vector2::zeros();

    let toi_separating = query::cast_shapes(
        &ball1_pos_intersecting,
        &ball1_vel_separating,
        &ball1,
        &ball2_pos,
        &ball2_vel,
        &ball2,
        ShapeCastOptions::default(),
    )
    .unwrap();

    let toi_penetrating_ignore_pen = query::cast_shapes(
        &ball1_pos_intersecting,
        &ball1_vel_penetrating,
        &ball1,
        &ball2_pos,
        &ball2_vel,
        &ball2,
        ShapeCastOptions {
            stop_at_penetration: false,
            ..Default::default()
        },
    )
    .unwrap();

    let toi_separating_ignore_pen = query::cast_shapes(
        &ball1_pos_intersecting,
        &ball1_vel_separating,
        &ball1,
        &ball2_pos,
        &ball2_vel,
        &ball2,
        ShapeCastOptions {
            stop_at_penetration: false,
            ..Default::default()
        },
    )
    .unwrap();

    assert_eq!(toi_separating.map(|hit| hit.time_of_impact), Some(0.0));

    assert_eq!(
        toi_penetrating_ignore_pen.map(|hit| hit.time_of_impact),
        Some(0.0),
    );

    assert_eq!(
        toi_separating_ignore_pen.map(|hit| hit.time_of_impact),
        None
    );
}

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

    let toi_intersecting = query::cast_shapes(
        &ball_pos_intersecting,
        &ball_vel1,
        &ball,
        &cuboid_pos,
        &cuboid_vel1,
        &cuboid,
        ShapeCastOptions::default(),
    )
    .unwrap();
    let toi_will_touch = query::cast_shapes(
        &ball_pos_will_touch,
        &ball_vel2,
        &ball,
        &cuboid_pos,
        &cuboid_vel2,
        &cuboid,
        ShapeCastOptions::default(),
    )
    .unwrap();
    let toi_wont_touch = query::cast_shapes(
        &ball_pos_wont_touch,
        &ball_vel2,
        &ball,
        &cuboid_pos,
        &cuboid_vel1,
        &cuboid,
        ShapeCastOptions::default(),
    )
    .unwrap();

    assert_eq!(toi_intersecting.map(|hit| hit.time_of_impact), Some(0.0));
    assert!(relative_eq!(
        toi_will_touch.unwrap().time_of_impact,
        ((2.0 as Real).sqrt() - 1.0) / (ball_vel2 - cuboid_vel2).norm()
    ));
    assert_eq!(toi_wont_touch.map(|hit| hit.time_of_impact), None);
}

#[test]
fn cuboid_cuboid_toi_issue_214() {
    let shape1 = Cuboid::new(Vector2::new(1.0, 1.0));
    let shape2 = Cuboid::new(Vector2::new(1.0, 1.5));

    let pos1 = Isometry2::new(Vector2::new(0.0, 0.0), 0.0);
    let pos2 = Isometry2::new(Vector2::new(10.0, 0.0), 0.0);

    let vel1 = Vector2::new(1.0, 0.0);
    let vel2 = Vector2::new(0.0, 0.0);

    let hit = query::cast_shapes(
        &pos1,
        &vel1,
        &shape1,
        &pos2,
        &vel2,
        &shape2,
        ShapeCastOptions::default(),
    )
    .unwrap();
    assert!(hit.is_some());
}

#[test]
fn cast_shapes_should_return_toi_for_ball_and_rotated_polyline() {
    let ball_isometry = Isometry2::identity();
    let ball_velocity = Vector2::new(1.0, 0.0);
    let ball = Ball::new(0.5);
    let polyline_isometry = Isometry2::rotation(-std::f32::consts::FRAC_PI_2);
    let polyline_velocity = Vector2::zeros();
    let polyline = Polyline::new(vec![Point2::new(1.0, 1.0), Point2::new(-1.0, 1.0)], None);

    assert_eq!(
        polyline_isometry.transform_point(&Point2::new(1.0, 1.0)),
        Point2::new(0.99999994, -1.0)
    );
    assert_eq!(
        polyline_isometry.transform_point(&Point2::new(-1.0, 1.0)),
        Point2::new(1.0, 0.99999994)
    );

    let hit = query::cast_shapes(
        &ball_isometry,
        &ball_velocity,
        &ball,
        &polyline_isometry,
        &polyline_velocity,
        &polyline,
        ShapeCastOptions::with_max_time_of_impact(1.0),
    )
    .unwrap();

    assert_eq!(hit.unwrap().time_of_impact, 0.5);
}

#[test]
fn cast_shapes_should_return_toi_for_ball_and_rotated_segment() {
    let ball_isometry = Isometry2::identity();
    let ball_velocity = Vector2::new(1.0, 0.0);
    let ball = Ball::new(0.5);
    let segment_isometry = Isometry2::rotation(-std::f32::consts::FRAC_PI_2);
    let segment_velocity = Vector2::zeros();
    let segment = Segment::new(Point2::new(1.0, 1.0), Point2::new(-1.0, 1.0));

    assert_eq!(
        segment_isometry.transform_point(&Point2::new(1.0, 1.0)),
        Point2::new(0.99999994, -1.0)
    );
    assert_eq!(
        segment_isometry.transform_point(&Point2::new(-1.0, 1.0)),
        Point2::new(1.0, 0.99999994)
    );

    let hit = query::cast_shapes(
        &ball_isometry,
        &ball_velocity,
        &ball,
        &segment_isometry,
        &segment_velocity,
        &segment,
        ShapeCastOptions::with_max_time_of_impact(1.0),
    )
    .unwrap();

    assert_eq!(hit.unwrap().time_of_impact, 0.49999994);
}

#[test]
fn cast_shapes_should_return_toi_for_rotated_segment_and_ball() {
    let ball_isometry = Isometry2::identity();
    let ball_velocity = Vector2::new(1.0, 0.0);
    let ball = Ball::new(0.5);
    let segment_isometry = Isometry2::rotation(-std::f32::consts::FRAC_PI_2);
    let segment_velocity = Vector2::zeros();
    let segment = Segment::new(Point2::new(1.0, 1.0), Point2::new(-1.0, 1.0));

    assert_eq!(
        segment_isometry.transform_point(&Point2::new(1.0, 1.0)),
        Point2::new(0.99999994, -1.0)
    );
    assert_eq!(
        segment_isometry.transform_point(&Point2::new(-1.0, 1.0)),
        Point2::new(1.0, 0.99999994)
    );

    let hit = query::cast_shapes(
        &segment_isometry,
        &segment_velocity,
        &segment,
        &ball_isometry,
        &ball_velocity,
        &ball,
        ShapeCastOptions::with_max_time_of_impact(1.0),
    )
    .unwrap();

    assert_eq!(hit.unwrap().time_of_impact, 0.5);
}
