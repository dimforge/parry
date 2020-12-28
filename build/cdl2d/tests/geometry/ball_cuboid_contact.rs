use cdl2d::query;
use cdl2d::shape::{Ball, Cuboid};
use nalgebra::{Isometry2, Vector2};
#[cfg(feature = "improved_fixed_point_support")]
use simba::scalar::FixedI40F24;

#[test]
fn test_ball_cuboid_query_contact() {
    let cuboid = Cuboid::new(Vector2::new(0.5, 0.5));
    let cuboid_pos = Isometry2::translation(0.0, 4.0);
    let ball = Ball::new(0.5);
    let ball_pos = Isometry2::translation(0.0517938, 3.05178815);
    let ct = query::contact(&cuboid_pos.inv_mul(&ball_pos), &cuboid, &ball, 0.0).unwrap();
    assert!(ct.is_some());
    let ct = query::contact(&ball_pos.inv_mul(&cuboid_pos), &ball, &cuboid, 0.0).unwrap();
    assert!(ct.is_some());
}
