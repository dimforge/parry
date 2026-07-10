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

#[test]
fn shape_cast_toi_accuracy_does_not_scale_with_shape_extent() {
    const BALL_RADIUS: Real = 0.5166;
    const MAX_TOI_ERROR: Real = 2.0e-4; // 0.2 mm.
    const MAX_WITNESS_ERROR: Real = 1.0e-6;
    const MAX_NORMAL_ERROR: Real = 5.0e-4;

    let ball = Ball::new(BALL_RADIUS);
    let direction = -Vector::Y;
    let mut toi_errors = Vec::new();

    for half_extent in [5.0, 50.0, 500.0] {
        let ground = Cuboid::new(Vector::new(half_extent, 0.5, half_extent));
        let ground_pose = Pose::translation(0.0, -0.5, 0.0);
        let mut max_toi_error: Real = 0.0;
        let mut max_witness_error: Real = 0.0;
        let mut max_normal_error: Real = 0.0;

        // Sweep sub-millimetre start offsets at several lateral positions. This mirrors
        // suspension probes near their resting pose while exercising large support points.
        for i in 0..200 {
            let y = 1.177 + i as Real * 5.0e-5;
            for (x, z) in [(1.4, 2.8), (-1.4, -0.4), (1.4, -2.0), (-1.4, 2.0)] {
                let ball_pose = Pose::translation(x, y, z);
                let hit = query::cast_shapes(
                    &ground_pose,
                    Vector::ZERO,
                    &ground,
                    &ball_pose,
                    direction,
                    &ball,
                    ShapeCastOptions::with_max_time_of_impact(2.0),
                )
                .unwrap()
                .expect("the downward cast should hit the cuboid");

                let expected_toi = y - BALL_RADIUS;
                max_toi_error = max_toi_error.max((hit.time_of_impact - expected_toi).abs());

                let witness1 = ground_pose.transform_point(hit.witness1);
                max_witness_error = max_witness_error.max(witness1.y.abs());
                max_normal_error = max_normal_error.max((hit.normal1 - Vector::Y).length());
            }
        }

        assert!(
            max_witness_error <= MAX_WITNESS_ERROR,
            "shape-cast witness error {max_witness_error} exceeded {MAX_WITNESS_ERROR} for half-extent {half_extent}",
        );
        assert!(
            max_normal_error <= MAX_NORMAL_ERROR,
            "shape-cast normal error {max_normal_error} exceeded {MAX_NORMAL_ERROR} for half-extent {half_extent}",
        );
        println!(
            "half-extent {half_extent}: max TOI error = {} mm, max witness error = {} mm",
            max_toi_error * 1000.0,
            max_witness_error * 1000.0,
        );
        toi_errors.push((half_extent, max_toi_error));
    }

    for (half_extent, max_toi_error) in toi_errors {
        assert!(
            max_toi_error <= MAX_TOI_ERROR,
            "shape-cast TOI error {max_toi_error} exceeded {MAX_TOI_ERROR} for half-extent {half_extent}",
        );
    }
}
