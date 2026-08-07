// Regression test for https://github.com/dimforge/parry/issues/373
//
// `cast_shapes_voxels_shape` delegates the cast to a cuboid centered on each candidate
// voxel, but stored `witness1` without shifting it back into the voxels shape's local
// frame, so it was off by the voxel's center (and `witness2` in the swapped wrapper).

use parry3d::math::{IVector, Pose, Vector};
use parry3d::query::{self, ShapeCastOptions};
use parry3d::shape::{Ball, Voxels};

fn floor_voxels() -> Voxels {
    // A single row of 20 unit voxels along +x: voxel `i` spans
    // [i, i + 1] x [0, 1] x [0, 1] in the shape's local space.
    let coords: Vec<_> = (0..20).map(|i| IVector::new(i, 0, 0)).collect();
    Voxels::new(Vector::new(1.0, 1.0, 1.0), &coords)
}

#[test]
fn voxels_shape_cast_witness1_is_in_voxels_local_space() {
    let voxels = floor_voxels();
    let ball = Ball::new(0.5);

    // Drop the ball straight down onto voxel 17, far from the shape's origin.
    let pos_voxels = Pose::identity();
    let pos_ball = Pose::translation(17.5, 3.5, 0.5);
    let vel_ball = Vector::new(0.0, -1.0, 0.0);

    let hit = query::cast_shapes(
        &pos_voxels,
        Vector::ZERO,
        &voxels,
        &pos_ball,
        vel_ball,
        &ball,
        ShapeCastOptions::default(),
    )
    .unwrap()
    .expect("the ball should hit the voxels floor");

    // Ball bottom starts at y = 3.0 and the voxel top face is at y = 1.0.
    assert!((hit.time_of_impact - 2.0).abs() < 1.0e-4);

    // `witness1` must lie on the hit voxel's top face, in the voxels shape's
    // local space (the bug reported it in the individual voxel's frame,
    // i.e. (0.0, 0.5, 0.0) here).
    let expected_witness1 = Vector::new(17.5, 1.0, 0.5);
    assert!(
        (hit.witness1 - expected_witness1).length() < 1.0e-2,
        "unexpected witness1: {:?}",
        hit.witness1
    );

    // The other fields were already correct: `witness2` is in the ball's frame
    // and the normals are translation-invariant.
    assert!((hit.witness2 - Vector::new(0.0, -0.5, 0.0)).length() < 1.0e-2);
    assert!((hit.normal1 - Vector::new(0.0, 1.0, 0.0)).length() < 1.0e-2);
    assert!((hit.normal2 - Vector::new(0.0, -1.0, 0.0)).length() < 1.0e-2);

    // The swapped shape order must agree with the direct one.
    let swapped_hit = query::cast_shapes(
        &pos_ball,
        vel_ball,
        &ball,
        &pos_voxels,
        Vector::ZERO,
        &voxels,
        ShapeCastOptions::default(),
    )
    .unwrap()
    .expect("the ball should hit the voxels floor");

    assert!((swapped_hit.time_of_impact - hit.time_of_impact).abs() < 1.0e-2);
    assert!((swapped_hit.witness1 - hit.witness2).length() < 1.0e-2);
    assert!((swapped_hit.witness2 - hit.witness1).length() < 1.0e-2);
    assert!((swapped_hit.normal1 - hit.normal2).length() < 1.0e-2);
    assert!((swapped_hit.normal2 - hit.normal1).length() < 1.0e-2);
}
