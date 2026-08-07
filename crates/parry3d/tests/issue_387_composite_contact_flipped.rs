// Regression test for https://github.com/dimforge/parry/issues/387
// (fixed by https://github.com/dimforge/parry/pull/388)
//
// `contact_composite_shape_shape` used to invert `pose12` a second time before
// traversing the composite's BVH, mirroring the traversal AABB and missing contacts for
// off-center parts. This checks it agrees with the flipped query.

use parry3d::math::{Pose, Vector};
use parry3d::query::contact;
use parry3d::shape::{Ball, Compound, SharedShape};

#[test]
fn contact_with_off_center_compound_part_agrees_with_flipped_order() {
    // A compound whose parts are far from its origin: a unit-half-extent cube
    // at (5, 0, 0) and another at (-5, -2, 0).
    let compound = Compound::new(vec![
        (
            Pose::translation(5.0, 0.0, 0.0),
            SharedShape::cuboid(1.0, 1.0, 1.0),
        ),
        (
            Pose::translation(-5.0, -2.0, 0.0),
            SharedShape::cuboid(1.0, 1.0, 1.0),
        ),
    ]);
    let ball = Ball::new(1.0);

    // The ball hovers 1.0 above the top face of the part at (5, 0, 0).
    let pos_compound = Pose::identity();
    let pos_ball = Pose::translation(5.0, 3.0, 0.0);
    let prediction = 1.5;

    let c12 = contact(&pos_compound, &compound, &pos_ball, &ball, prediction)
        .unwrap()
        .expect("the ball is within the prediction distance of the off-center part");

    // With the double inversion, the traversal AABB was around (-5, -3, 0):
    // no part there, hence no contact (or a contact against the wrong part).
    assert!((c12.dist - 1.0).abs() < 1.0e-3, "dist = {}", c12.dist);
    assert!((c12.point1 - Vector::new(5.0, 1.0, 0.0)).length() < 1.0e-3);
    assert!((c12.point2 - Vector::new(5.0, 2.0, 0.0)).length() < 1.0e-3);
    assert!((c12.normal1 - Vector::new(0.0, 1.0, 0.0)).length() < 1.0e-3);
    assert!((c12.normal2 - Vector::new(0.0, -1.0, 0.0)).length() < 1.0e-3);

    // The flipped order must return the same contact with the roles swapped.
    let c21 = contact(&pos_ball, &ball, &pos_compound, &compound, prediction)
        .unwrap()
        .expect("the flipped query must detect the same contact");

    assert!((c21.dist - c12.dist).abs() < 1.0e-6);
    assert!((c21.point1 - c12.point2).length() < 1.0e-6);
    assert!((c21.point2 - c12.point1).length() < 1.0e-6);
    assert!((c21.normal1 - c12.normal2).length() < 1.0e-6);
    assert!((c21.normal2 - c12.normal1).length() < 1.0e-6);
}
