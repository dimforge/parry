use cdl3d::query::{time_of_impact, DefaultQueryDispatcher};
use cdl3d::shape::Cuboid;
use na::{self, Isometry3, Vector3};

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
    let pos1 = Isometry3::translation(0.0, 1.1, 0.0);
    let pos2 = Isometry3::identity();
    let vel1 = Vector3::y() * v_y;
    let vel2 = Vector3::zeros();
    let cuboid = Cuboid::new(Vector3::new(0.5, 0.5, 0.5));

    time_of_impact(
        &(pos1.inverse() * pos2),
        &(vel2 - vel1),
        &cuboid,
        &cuboid,
        std::f32::MAX,
        0.0,
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
