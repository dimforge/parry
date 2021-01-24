use na::{self, Isometry2, Vector2};
use parry2d::query::{self, ContactManifold, DefaultQueryDispatcher, PersistentQueryDispatcher};
use parry2d::shape::Cuboid;

#[test]
#[allow(non_snake_case)]
fn cuboid_cuboid_EPA() {
    let c = Cuboid::new(Vector2::new(2.0, 1.0));
    let m1 = Isometry2::translation(3.5, 0.0);
    let m2 = Isometry2::identity();

    let res = query::details::contact_support_map_support_map(&m1.inv_mul(&m2), &c, &c, 10.0)
        .expect("Penetration not found.");
    assert_eq!(res.dist, -0.5);
    assert_eq!(res.normal1, -Vector2::x_axis());

    let m1 = Isometry2::translation(0.0, 0.2);
    let res = query::details::contact_support_map_support_map(&m1.inv_mul(&m2), &c, &c, 10.0)
        .expect("Penetration not found.");
    assert_eq!(res.dist, -1.8);
    assert_eq!(res.normal1, -Vector2::y_axis());
}

#[test]
fn cuboids_large_size_ratio_issue_181() {
    let cuboid_a = Cuboid::new(Vector2::new(10.0, 10.0));
    let cuboid_b = Cuboid::new(Vector2::new(300.0, 1.5));

    let pos_b = Isometry2::new(Vector2::new(5.0, 0.0), 1.5);

    let dispatcher = DefaultQueryDispatcher;
    let mut p = Vector2::new(0.0, 0.0);
    let mut angle = 0.0;

    // Used to panic at some point:
    // thread 'main' panicked at 'assertion failed: neg_dist <= gjk::eps_tol()', parry_geometry/query/algorithms/EPA.rs:26:9
    for _ in 1..200000 {
        p.x += 0.0001;
        angle += 0.005;

        let pos_a = Isometry2::new(p, angle);
        let pos_ab = pos_a.inv_mul(&pos_b);
        let mut manifold: ContactManifold<(), ()> = ContactManifold::new();
        dispatcher
            .contact_manifold_convex_convex(&pos_ab, &cuboid_a, &cuboid_b, 0.0, &mut manifold)
            .unwrap();

        if let Some(deepest) = manifold.find_deepest_contact() {
            p += pos_a * manifold.local_n1 * deepest.dist;
        }
    }
}
