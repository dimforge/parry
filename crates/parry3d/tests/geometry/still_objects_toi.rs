use parry3d::math::{GlamVectorOps, Isometry, Vector};
use parry3d::query::time_of_impact;
use parry3d::shape::Cuboid;

/**
 * Issue #141
 * Sets up a situation like this:
 * ```raw
 * +---+
 * | 1 |
 * |   |
 * +---+
 *
 * +---+
 * | 2 |
 * |   |
 * +---+
 * ```
 * with box 1 having the provided v_y.
 */
fn collide(v_y: f32) -> Option<f32> {
    let pos1 = Isometry::translation(0.0, 1.1, 0.0);
    let pos2 = Isometry::identity();
    let vel1 = Vector::y() * v_y;
    let vel2 = Vector::zeros();
    let cuboid = Cuboid::new(Vector::new(0.5, 0.5, 0.5));

    time_of_impact(
        &pos1,
        &vel1,
        &cuboid,
        &pos2,
        &vel2,
        &cuboid,
        std::f32::MAX,
        true,
    )
    .unwrap()
    .map(|toi| toi.toi)
}

#[test]
fn no_movement() {
    assert_eq!(collide(0.0), None);
}

#[test]
fn moving_up_misses() {
    assert_eq!(collide(1.0), None);
}

#[test]
fn moving_down_hits() {
    assert!(collide(-1.0).is_some());
}
