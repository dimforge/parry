//! Regression tests for shape-cast TOI against very large cuboids.
//!
//! Before the fix in `minkowski_ray_cast`, sweeping a cylinder straight
//! down onto a 10 000 × 1 × 10 000 "backstop" cuboid returned ~10× too
//! small a `time_of_impact`.

use parry3d::math::{Pose, Real, Vector};
use parry3d::query::{self, ShapeCastOptions};
use parry3d::shape::{Ball, Cuboid, Cylinder};

fn assert_close(actual: Real, expected: Real, tol: Real, label: &str) {
    assert!(
        (actual - expected).abs() <= tol,
        "{}: got {}, expected {} (tol {})",
        label,
        actual,
        expected,
        tol,
    );
}

fn down_cast_options() -> ShapeCastOptions {
    ShapeCastOptions {
        max_time_of_impact: 2.0,
        target_distance: 0.0,
        stop_at_penetration: true,
        compute_impact_geometry_on_penetration: true,
    }
}

#[test]
fn ball_small_cuboid_down_sweep() {
    let b = Ball::new(0.438_f32);
    let cub = Cuboid::new(Vector::new(1.0, 0.5, 1.0));

    let hit = query::cast_shapes(
        &Pose::translation(0.0, 1.281, 0.0),
        Vector::new(0.0, -1.0, 0.0),
        &b,
        &Pose::translation(0.0, -0.6, 0.0),
        Vector::ZERO,
        &cub,
        down_cast_options(),
    )
    .unwrap()
    .expect("should hit");

    assert_close(hit.time_of_impact, 0.943, 1e-3, "ball_small_cuboid TOI");
}

#[test]
fn cylinder_small_cuboid_down_sweep() {
    let cyl = Cylinder::new(0.25_f32, 0.438);
    let cub = Cuboid::new(Vector::new(1.0, 0.5, 1.0));

    // Bottom cap y = 1.281 - 0.25 = 1.031; cuboid top y = -0.1 → TOI = 1.131.
    let hit = query::cast_shapes(
        &Pose::translation(0.0, 1.281, 0.0),
        Vector::new(0.0, -1.0, 0.0),
        &cyl,
        &Pose::translation(0.0, -0.6, 0.0),
        Vector::ZERO,
        &cub,
        down_cast_options(),
    )
    .unwrap()
    .expect("should hit");

    assert_close(hit.time_of_impact, 1.131, 1e-3, "cyl_small_cuboid TOI");
}

#[test]
fn ball_huge_cuboid_down_sweep() {
    let b = Ball::new(0.438_f32);
    let cub = Cuboid::new(Vector::new(5000.0, 0.5, 5000.0));

    let hit = query::cast_shapes(
        &Pose::translation(43.0, 1.281, 5.7),
        Vector::new(0.0, -1.0, 0.0),
        &b,
        &Pose::translation(0.0, -0.6, 0.0),
        Vector::ZERO,
        &cub,
        down_cast_options(),
    )
    .unwrap()
    .expect("should hit");

    assert_close(hit.time_of_impact, 0.943, 1e-3, "ball_huge_cuboid TOI");
}

// Pre-fix: `time_of_impact` ≈ 0.12 instead of ~1.13.
#[test]
fn cylinder_huge_cuboid_down_sweep() {
    let cyl = Cylinder::new(0.25_f32, 0.438);
    let cub = Cuboid::new(Vector::new(5000.0, 0.5, 5000.0));

    let hit = query::cast_shapes(
        &Pose::translation(43.0, 1.281, 5.7),
        Vector::new(0.0, -1.0, 0.0),
        &cyl,
        &Pose::translation(0.0, -0.6, 0.0),
        Vector::ZERO,
        &cub,
        down_cast_options(),
    )
    .unwrap()
    .expect("should hit");

    // Loose tolerance: GJK precision is inherently limited on f32 supports
    // with ±5000 magnitude. The pre-fix result was off by ~10×.
    assert_close(hit.time_of_impact, 1.131, 5e-3, "cyl_huge_cuboid TOI");
}
