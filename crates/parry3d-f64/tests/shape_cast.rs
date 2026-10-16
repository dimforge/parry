use parry3d_f64::math::{Pose, Vector};
use parry3d_f64::query::{cast_shapes, ShapeCastOptions};
use parry3d_f64::shape::Cuboid;

#[test]
fn cast_starting_at_touching_boundary_and_moving_inward_hits() {
    let floor = Cuboid::new(Vector::new(38.0, 0.5, 4.0));
    let floor_pos = Pose::from_translation(Vector::new(42.0, -0.5, 2.0));
    let boxy = Cuboid::new(Vector::splat(0.25));
    let options = ShapeCastOptions {
        max_time_of_impact: 6.54806820868758261e-3,
        target_distance: 0.01,
        stop_at_penetration: false,
        compute_impact_geometry_on_penetration: true,
    };
    let start = Vector::new(
        f64::from_bits(0x4038019999998b41),
        f64::from_bits(0x3fd0a3d70a3d7099),
        f64::from_bits(0xbe095faa25f26d8d),
    );
    let dir = Vector::new(9.54479978035029708e-1, -2.98274993135946798e-1, 0.0);

    for dy in [0.0, -1e-9, 1e-9] {
        let pos = Pose::from_translation(start + Vector::Y * dy);
        let hit = cast_shapes(&floor_pos, Vector::ZERO, &floor, &pos, dir, &boxy, options).unwrap();
        assert!(hit.is_some(), "cast missed with start offset {dy}");
    }
}
