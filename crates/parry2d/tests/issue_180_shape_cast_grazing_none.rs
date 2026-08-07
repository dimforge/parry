// Regression test for https://github.com/dimforge/parry/issues/180
// (MRE from the issue + the capsule-vs-segment repro of PR #410)
//
// Grazing casts returned `None` on 1-ulp-sensitive inputs: the GJK Minkowski ray-cast
// reached a full-dimensional simplex with `min_bound` barely above the absolute
// tolerance and rejected the valid hit found so far.

use parry2d::math::{Pose, Vector};
use parry2d::query::details::cast_shapes_support_map_support_map;
use parry2d::query::{self, ShapeCastOptions};
use parry2d::shape::{Capsule, Cuboid, Segment};

fn cast(pos12: Pose, vel12: Vector, segment: &Segment, cuboid: &Cuboid) -> bool {
    cast_shapes_support_map_support_map(&pos12, vel12, segment, cuboid, ShapeCastOptions::default())
        .is_some()
}

#[test]
fn grazing_cast_is_not_ulp_sensitive() {
    let vel = Vector::new(0.0081385, -0.999966);
    let segment = Segment::new(Vector::new(1.0, 1.0), Vector::new(124.0, 1.0));
    let cuboid = Cuboid::new(Vector::new(2.0, 2.0));

    // These two always worked.
    assert!(cast(
        Pose::translation(15.689464, 54.00709),
        vel,
        &segment,
        &cuboid
    ));
    assert!(cast(
        Pose::translation(15.689466, 54.00709),
        vel,
        &segment,
        &cuboid
    ));

    // This one, 1 ulp away from the previous two, used to return `None`.
    assert!(cast(
        Pose::translation(15.689465, 54.00709),
        vel,
        &segment,
        &cuboid
    ));
}

// From PR #410: a capsule 0.1 units away from a steep segment, moving straight
// toward it, must produce a hit.
#[test]
fn capsule_segment_steep_slope_toi() {
    // Segment in local space (centered at entity transform [-648, -288]).
    let segment = Segment::new(Vector::new(-24.0, 48.0), Vector::new(24.0, -48.0));
    let capsule = Capsule::new_y(14.0, 8.0);

    // Exact game positions.
    let segment_pose = Pose::translation(-648.0, -288.0);
    let capsule_pose = Pose::translation(-653.3891, -245.08746);

    // Per-frame remaining velocity (vel / 60fps).
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

    assert!(
        result.is_some(),
        "GJK missed a collision: capsule is 0.1 units from a steep segment \
         and moving directly toward it. This is a clear hit that should \
         never be missed.",
    );
}
