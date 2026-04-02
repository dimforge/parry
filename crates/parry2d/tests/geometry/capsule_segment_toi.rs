use parry2d::math::{Pose, Vector};
use parry2d::query::{self, ShapeCastOptions};
use parry2d::shape::{Capsule, Segment};

/// A capsule 0.1 units from a steep segment should detect collision when
/// cast toward it. Uses exact coordinates from a game scenario where
/// GJK's minkowski_ray_cast incorrectly returns None.
///
/// The segment [-672,-240] -> [-624,-336] is a 63.4 steep wall.
/// The capsule (radius=8, half_height=14) at [-653.39, -245.09] is
/// 0.1 units from the surface and moving toward it at a velocity of [-200, -325]/frame.
///
/// Bug: minkowski_ray_cast reaches a full simplex with min_bound slightly
/// above eps_tol while ltoi > 0, and returns None despite having found
/// valid contact.
#[test]
fn capsule_segment_steep_slope_toi() {
    // Segment in local space (centered at entity transform [-648, -288])
    let segment = Segment::new(Vector::new(-24.0, 48.0), Vector::new(24.0, -48.0));
    let capsule = Capsule::new_y(14.0, 8.0);

    // Exact game positions
    let segment_pose = Pose::new(Vector::new(-648.0, -288.0), 0.0);
    let capsule_pose = Pose::new(Vector::new(-653.3891, -245.08746), 0.0);

    // Per-frame remaining velocity (vel / 60fps)
    let capsule_vel = Vector::new(-3.3333337, -5.4166675);

    let options = ShapeCastOptions {
        max_time_of_impact: capsule_vel.length(),
        ..Default::default()
    };

    let result = query::cast_shapes(
        &capsule_pose,
        capsule_vel,
        &capsule,
        &segment_pose,
        Vector::ZERO,
        &segment,
        options,
    )
    .unwrap();

    // The capsule is 0.1 units from the segment surface, moving toward it.
    // The collision should be detected well within the cast distance.
    assert!(
        result.is_some(),
        "GJK missed a collision: capsule is 0.1 units from a steep segment \
         and moving directly toward it. This is a clear hit that should \
         never be missed.",
    );
}
