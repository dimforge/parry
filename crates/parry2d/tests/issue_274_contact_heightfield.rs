// Regression test for https://github.com/dimforge/parry/issues/274 (2D variant)
//
// `query::contact` with a heightfield used to return `Err(Unsupported)`.

use parry2d::math::{Pose, Vector};
use parry2d::query::contact;
use parry2d::shape::{Ball, HeightField};

const EPS: f32 = 1.0e-4;

#[test]
fn ball_vs_heightfield_both_orders() {
    // A flat heightfield at y = 0 spanning x in [-5, 5].
    let hf = HeightField::new(vec![0.0, 0.0, 0.0], Vector::new(10.0, 1.0));
    let ball = Ball::new(0.5);
    let pos_hf = Pose::identity();
    let pos_ball = Pose::from_translation(Vector::new(1.0, 0.3));

    let c = contact(&pos_hf, &hf, &pos_ball, &ball, 0.0)
        .unwrap()
        .unwrap();
    assert!(
        (c.dist + 0.2).abs() < EPS,
        "expected dist ~ -0.2, got {}",
        c.dist
    );
    assert!(
        (c.normal1 - Vector::Y).length() < EPS,
        "normal1 should be +Y, got {:?}",
        c.normal1
    );
    assert!((c.point1 - Vector::new(1.0, 0.0)).length() < EPS);

    let c = contact(&pos_ball, &ball, &pos_hf, &hf, 0.0)
        .unwrap()
        .unwrap();
    assert!((c.dist + 0.2).abs() < EPS);
    assert!(
        (c.normal1 + Vector::Y).length() < EPS,
        "normal1 should be -Y, got {:?}",
        c.normal1
    );

    // Separated within prediction.
    let pos_ball = Pose::from_translation(Vector::new(-2.0, 1.0));
    let c = contact(&pos_hf, &hf, &pos_ball, &ball, 1.0)
        .unwrap()
        .unwrap();
    assert!(
        (c.dist - 0.5).abs() < EPS,
        "expected dist ~ 0.5, got {}",
        c.dist
    );

    // Separated beyond prediction: None, not Unsupported.
    let c = contact(&pos_hf, &hf, &pos_ball, &ball, 0.1).unwrap();
    assert!(c.is_none());
}
