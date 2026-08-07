// Regression test for https://github.com/dimforge/parry/issues/396
//
// `intersection_test` returned `Ok(false)` for two overlapping coaxial cylinders: their
// degenerate axial support directions stall the GJK simplex, and the stagnation exits
// reported `Proximity` (=> disjoint) without ever certifying separation.

use parry3d::math::{Pose, Vector};
use parry3d::query::intersection_test;
use parry3d::shape::{Capsule, Cylinder};

#[test]
fn coaxial_overlapping_cylinders_intersect() {
    // The exact setup from the issue.
    let iso = Pose::IDENTITY;
    let big = Cylinder::new(50.0, 100.0);
    let small = Cylinder::new(1.5, 1.0);

    assert_eq!(intersection_test(&iso, &big, &iso, &small), Ok(true));
}

#[test]
fn coaxial_cylinder_capsule_intersect() {
    // The second setup from the issue (cylinder vs tall thin capsule).
    let iso = Pose::IDENTITY;
    let big = Cylinder::new(50.0, 100.0);
    let small = Capsule::new(
        Vector::new(0.0, -50.0, 0.0),
        Vector::new(0.0, 50.0, 0.0),
        10.0,
    );

    assert_eq!(intersection_test(&iso, &big, &iso, &small), Ok(true));
}

#[test]
fn offset_and_rotated_cylinders_intersecting() {
    let big = Cylinder::new(50.0, 100.0);
    let small = Cylinder::new(1.5, 1.0);

    // Slightly offset (but still deeply overlapping) pairs.
    for offset in [
        Vector::new(0.01, 0.0, 0.0),
        Vector::new(0.0, 0.01, 0.0),
        Vector::new(0.0, 0.0, -0.01),
        Vector::new(3.0, 5.0, -2.0),
    ] {
        let pos2 = Pose::from_translation(offset);
        assert_eq!(
            intersection_test(&Pose::IDENTITY, &big, &pos2, &small),
            Ok(true),
            "offset {offset:?} should intersect"
        );
    }

    // Rotated overlapping pairs.
    for angle in [0.01, 0.5, core::f32::consts::FRAC_PI_2] {
        let pos2 = Pose::rotation(Vector::new(0.0, 0.0, angle));
        assert_eq!(
            intersection_test(&Pose::IDENTITY, &big, &pos2, &small),
            Ok(true),
            "rotation angle {angle} should intersect"
        );
    }
}

#[test]
fn coaxial_cylinders_separated() {
    // Guard the other direction: separated (still coaxial/axisymmetric) pairs
    // must keep returning `false`.
    let big = Cylinder::new(50.0, 100.0);
    let small = Cylinder::new(1.5, 1.0);

    // Separated along the symmetry axis.
    for dy in [51.6, 60.0, 200.0] {
        let pos2 = Pose::translation(0.0, dy, 0.0);
        assert_eq!(
            intersection_test(&Pose::IDENTITY, &big, &pos2, &small),
            Ok(false),
            "axial offset {dy} should be disjoint"
        );
    }

    // Separated radially.
    for dx in [101.1, 110.0, 500.0] {
        let pos2 = Pose::translation(dx, 0.0, 0.0);
        assert_eq!(
            intersection_test(&Pose::IDENTITY, &big, &pos2, &small),
            Ok(false),
            "radial offset {dx} should be disjoint"
        );
    }

    // Separated and rotated.
    let pos2 = Pose::new(Vector::new(0.0, 53.0, 0.0), Vector::new(0.0, 0.0, 0.7));
    assert_eq!(
        intersection_test(&Pose::IDENTITY, &big, &pos2, &small),
        Ok(false),
        "rotated separated pair should be disjoint"
    );
}
