// Regression test for https://github.com/dimforge/parry/issues/106
//
// A shape-cast starting exactly touching (`toi == 0`, not penetrating) returned
// `Converged` for ball-vs-ball but `PenetratingOrWithinTargetDist` for the generic
// support-map path, which flagged any zero TOI as penetrating.

use parry3d::math::{Pose, Vector};
use parry3d::query::{self, ShapeCastOptions, ShapeCastStatus};
use parry3d::shape::{Ball, Cuboid};

fn cast_status(
    pos1: Pose,
    g1: &impl parry3d::shape::Shape,
    pos2: Pose,
    g2: &impl parry3d::shape::Shape,
    options: ShapeCastOptions,
) -> ShapeCastStatus {
    query::cast_shapes(
        &pos1,
        Vector::new(0.0, -1.0, 0.0),
        g1,
        &pos2,
        Vector::ZERO,
        g2,
        options,
    )
    .unwrap()
    .expect("the cast should return a hit")
    .status
}

#[test]
fn touching_at_toi_zero_is_converged_for_all_backends() {
    let ball = Ball::new(0.5);
    let other_ball = Ball::new(0.5);
    let cuboid = Cuboid::new(Vector::new(1.0, 1.0, 1.0));

    // Exactly touching configurations (all coordinates exact in floats):
    // - ball at y = 1.0 touching a ball at the origin (r1 + r2 = 1.0);
    // - ball at y = 1.5 touching the top face (y = 1.0) of the cuboid.
    let ball_ball_status = cast_status(
        Pose::translation(0.0, 1.0, 0.0),
        &ball,
        Pose::IDENTITY,
        &other_ball,
        ShapeCastOptions::default(),
    );
    let ball_cuboid_status = cast_status(
        Pose::translation(0.0, 1.5, 0.0),
        &ball,
        Pose::IDENTITY,
        &cuboid,
        ShapeCastOptions::default(),
    );

    assert_eq!(ball_ball_status, ShapeCastStatus::Converged);
    assert_eq!(ball_cuboid_status, ball_ball_status);

    // Same, without the impact-geometry fallback.
    let no_geometry_options = ShapeCastOptions {
        compute_impact_geometry_on_penetration: false,
        ..Default::default()
    };
    let ball_ball_status = cast_status(
        Pose::translation(0.0, 1.0, 0.0),
        &ball,
        Pose::IDENTITY,
        &other_ball,
        no_geometry_options,
    );
    let ball_cuboid_status = cast_status(
        Pose::translation(0.0, 1.5, 0.0),
        &ball,
        Pose::IDENTITY,
        &cuboid,
        no_geometry_options,
    );
    assert_eq!(ball_ball_status, ShapeCastStatus::Converged);
    assert_eq!(ball_cuboid_status, ball_ball_status);
}

#[test]
fn penetrating_at_toi_zero_is_reported_for_all_backends() {
    let ball = Ball::new(0.5);
    let other_ball = Ball::new(0.5);
    let cuboid = Cuboid::new(Vector::new(1.0, 1.0, 1.0));

    // Actually overlapping by 0.125.
    let ball_ball_status = cast_status(
        Pose::translation(0.0, 0.875, 0.0),
        &ball,
        Pose::IDENTITY,
        &other_ball,
        ShapeCastOptions::default(),
    );
    let ball_cuboid_status = cast_status(
        Pose::translation(0.0, 1.375, 0.0),
        &ball,
        Pose::IDENTITY,
        &cuboid,
        ShapeCastOptions::default(),
    );

    assert_eq!(
        ball_ball_status,
        ShapeCastStatus::PenetratingOrWithinTargetDist
    );
    assert_eq!(ball_cuboid_status, ball_ball_status);
}
