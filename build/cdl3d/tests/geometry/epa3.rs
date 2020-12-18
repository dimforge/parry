use cdl3d::query;
use cdl3d::shape::Cuboid;
use na::{self, Isometry3, Vector3};

#[test]
#[allow(non_snake_case)]
fn cuboid_cuboid_EPA() {
    let c = Cuboid::new(Vector3::new(2.0, 1.0, 1.0));
    let m1 = Isometry3::translation(3.5, 0.0, 0.0);
    let m2 = Isometry3::identity();

    let res = query::details::contact_support_map_support_map(&m1.inv_mul(&m2), &c, &c, 10.0)
        .expect("Penetration not found.");
    assert_eq!(res.dist, -0.5);
    assert_eq!(res.normal1, -Vector3::x_axis());

    let m1 = Isometry3::translation(0.0, 0.2, 0.0);
    let res = query::details::contact_support_map_support_map(&m1.inv_mul(&m2), &c, &c, 10.0)
        .expect("Penetration not found.");
    assert_eq!(res.dist, -1.8);
    assert_eq!(res.normal1, -Vector3::y_axis());
}
