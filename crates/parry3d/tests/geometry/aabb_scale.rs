use parry3d::bounding_volume::Aabb;
use parry3d::math::Vector;

#[test]
fn test_aabb_scale_wrt_center() {
    let aabb = Aabb::from_half_extents(Vector::new(1.0, 2.0, 3.0), Vector::new(4.0, 5.0, 6.0));
    let scale = Vector::new(10.0, -20.0, 50.0);
    let scaled_aabb = aabb.scaled_wrt_center(scale);
    let scaled_aabb_neg = aabb.scaled_wrt_center(-scale);
    let scaled_aabb_abs = aabb.scaled_wrt_center(scale.abs());

    assert_eq!(&scaled_aabb, &scaled_aabb_neg);
    assert_eq!(&scaled_aabb, &scaled_aabb_abs);
    assert_eq!(aabb.center(), scaled_aabb.center());
    assert_eq!(scaled_aabb.half_extents(), Vector::new(40.0, 100.0, 300.0));
}
