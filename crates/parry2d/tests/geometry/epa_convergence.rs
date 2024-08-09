use na::Vector2;
use parry2d::{
    math::{Isometry, Point, Real},
    query,
    shape::{Capsule, ConvexPolygon, SharedShape},
};

/// Original issue: https://github.com/dimforge/parry/issues/205
#[test]
fn capsule_convergence() {
    let shape1 = Capsule::new_y(5.0, 10.0);
    let mut vec = Vec::<Point<Real>>::with_capacity(3);
    vec.push(Point::<Real>::new(64.0, 507.0));
    vec.push(Point::<Real>::new(440.0, 326.0));
    vec.push(Point::<Real>::new(1072.0, 507.0));
    let shape2 = ConvexPolygon::from_convex_polyline(vec);
    let shape2 = shape2.unwrap();
    let transform1 = Isometry::new(Vector2::new(381.592, 348.491), 0.0);
    let transform2 = Isometry::new(Vector2::new(0.0, 0.0), 0.0);

    let _res = query::details::contact_support_map_support_map(
        &transform1.inv_mul(&transform2),
        &shape1,
        &shape2,
        10.0,
    )
    .expect("Penetration not found.");
    let shared_shape1 = SharedShape::new(shape1);
    let shared_shape2 = SharedShape::new(shape2);

    if let Ok(Some(_contact)) = query::contact(
        &transform1,
        shared_shape1.as_ref(),
        &transform2,
        shared_shape2.as_ref(),
        1.0,
    ) {
        println!("collision");
    } else {
        panic!("no collision");
    }
}
