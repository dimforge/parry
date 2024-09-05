use na::{self, Isometry3, Point3, Vector3};
use parry3d::query;
use parry3d::query::gjk::VoronoiSimplex;
use parry3d::shape::{Cuboid, Triangle};

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

#[test]
fn triangle_vertex_touches_triangle_edge_epa() {
    // Related issues:
    // https://github.com/dimforge/parry/issues/253
    // https://github.com/dimforge/parry/issues/246

    let mesh1 = Triangle::new(
        Point3::new(-13.174434, 1.0, 8.736801),
        Point3::new(3.5251038, 1.0, 12.1),
        Point3::new(3.2048466, 1.0, 12.218325),
    );
    let mesh2 = Triangle::new(
        Point3::new(-1.63, 0.0, 11.19),
        Point3::new(-2.349647, 0.0, 11.037681),
        Point3::new(-2.349647, 1.0, 11.037681),
    );

    // TODO: Check the return-value of the function once
    // it the function's correctness issue have been taken care of
    // Right now we just want it not to crash
    let gjk_result = query::details::contact_support_map_support_map_with_params(
        &Isometry3::identity(),
        &mesh1,
        &mesh2,
        0.00999999977,
        &mut VoronoiSimplex::new(),
        None,
    );

    assert!(!matches!(gjk_result, query::gjk::GJKResult::NoIntersection(_)), "PARTIAL SUCCESS: contact_support_map_support_map_with_params did not crash but did not produce the desired result");
}