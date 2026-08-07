// Regression test for https://github.com/dimforge/rapier/issues/810
//
// A small cube resting face-down on a large thin cylinder disc could get a single-point
// manifold: the cap's circular face is approximated by an inscribed square anchored at
// an arbitrary azimuth, missing up to ~36% of the cap near its rim. The approximation is
// now oriented toward the contact point, so face contacts clip to a full manifold.

use parry3d::math::{Pose, Real, Vector};
use parry3d::query::{ContactManifold, DefaultQueryDispatcher, PersistentQueryDispatcher};
use parry3d::shape::{Cuboid, Cylinder};

const CAP_RADIUS: Real = 10.0;

fn cap_manifold(x: Real, z: Real, sep: Real) -> ContactManifold<(), ()> {
    // The issue's geometry: thin disc (radius 10, half-height 0.05), small cube
    // (half-extents 0.05) face-down on the top cap.
    let disc = Cylinder::new(0.05, CAP_RADIUS);
    let cube = Cuboid::new(Vector::splat(0.05));
    let pos12 = Pose::translation(x, 0.05 + 0.05 + sep, z);
    let mut manifold = ContactManifold::new();
    DefaultQueryDispatcher
        .contact_manifold_convex_convex(&pos12, &disc, &cube, None, None, 0.01, &mut manifold)
        .unwrap();
    manifold
}

/// A cube face resting anywhere on the cap must produce a multi-point manifold with a
/// near-axial normal — including in the band near the rim not covered by an arbitrarily
/// anchored inscribed square (`|x| + |z| > radius`, e.g. the issue's tunneling cubes).
#[test]
fn cuboid_on_cylinder_cap_has_multi_point_manifold() {
    for &(x, z) in &[
        (0.0, 0.0),                 // cap center
        (3.0, -2.0),                // mid cap
        (5.54, 4.51),               // rim band: |x| + |z| > radius (used to get 1 point)
        (-5.65, 5.65),              // worst-case azimuth on the rim band
        (7.0, 7.0),                 // radius ~9.9, close to the rim
        (9.9, 0.0),                 // right at the rim
    ] {
        for &sep in &[-0.001, 0.0001] {
            let manifold = cap_manifold(x, z, sep);
            assert!(
                manifold.points.len() >= 3,
                "cube at ({x}, {z}), sep {sep}: got {} contact point(s), need >= 3 \
                 for a stable face-on-cap contact",
                manifold.points.len()
            );
            assert!(
                manifold.local_n1.y.abs() > 0.99,
                "cube at ({x}, {z}), sep {sep}: expected a near-axial normal, got {:?}",
                manifold.local_n1
            );
        }
    }
}

/// A cube touching the cylinder's curved side must keep the current behavior: a radial
/// normal, with every contact point on the side's support segment (the line of the
/// cylinder's curved surface closest to the cube).
#[test]
fn cuboid_on_cylinder_side_unchanged() {
    let cylinder = Cylinder::new(2.0, 0.5);
    let cube = Cuboid::new(Vector::splat(0.25));
    let pos12 = Pose::translation(0.5 + 0.25 - 0.001, 0.0, 0.0);
    let mut manifold = ContactManifold::<(), ()>::new();
    DefaultQueryDispatcher
        .contact_manifold_convex_convex(&pos12, &cylinder, &cube, None, None, 0.01, &mut manifold)
        .unwrap();
    assert!(!manifold.points.is_empty());
    assert!(
        manifold.local_n1.x > 0.99,
        "expected a radial normal, got {:?}",
        manifold.local_n1
    );
    for pt in &manifold.points {
        assert!(
            (pt.local_p1.x - 0.5).abs() < 1.0e-3 && pt.local_p1.z.abs() < 1.0e-3,
            "side contact point should lie on the cylinder's support segment, got {:?}",
            pt.local_p1
        );
    }
}

/// A cylinder standing on its cap on a cuboid floor must also get a multi-point manifold.
#[test]
fn cylinder_standing_on_cuboid_floor() {
    let floor = Cuboid::new(Vector::new(5.0, 0.5, 5.0));
    let cylinder = Cylinder::new(0.4, 0.3);
    let pos12 = Pose::translation(1.0, 0.5 + 0.4 - 0.001, -2.0);
    let mut manifold = ContactManifold::<(), ()>::new();
    DefaultQueryDispatcher
        .contact_manifold_convex_convex(&pos12, &floor, &cylinder, None, None, 0.01, &mut manifold)
        .unwrap();
    assert!(
        manifold.points.len() >= 3,
        "standing cylinder should rest on >= 3 cap points, got {}",
        manifold.points.len()
    );
    assert!(manifold.local_n1.y > 0.99);
}
