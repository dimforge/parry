// Verification tests for https://github.com/dimforge/parry/issues/157
// (ported from Davidster's fork; frustum shapes pulled from a game engine)
//
// `contact`/`intersection_test` between two view-frustum polyhedra used to fail: GJK
// reported the intersection but EPA hit its iteration cap and returned `None`. EPA now
// caps at 100 iterations and returns its best-effort result.

use parry3d::math::{Pose, Vector};
use parry3d::query;
use parry3d::query::PointQuery;
use parry3d::shape::ConvexPolyhedron;

#[test]
fn convex_polyhedra_contact() {
    let convex_polyhedron_a_points = [
        [-0.7938391, 5.1101756, 0.1773476],
        [-5.293839, 0.6101756, -4.3226523],
        [-5.293839, 0.6101756, 5.6773477],
        [-0.7938391, 5.1101756, 1.1773477],
        [-0.7938391, 6.1101756, 0.1773476],
        [-5.293839, 10.610176, -4.3226523], // this point and only this point lies inside of b
        [-5.293839, 10.610176, 5.6773477],
        [-0.7938391, 6.1101756, 1.1773477],
    ]
    .map(|arr| Vector::new(arr[0], arr[1], arr[2]));
    let convex_polyhedron_a =
        ConvexPolyhedron::from_convex_hull(&convex_polyhedron_a_points).unwrap();

    let convex_polyhedron_b_points = [
        [8.114634, 4.6308937, 0.76987207],
        [-18.83664, -53.96347, -649.29565],
        [-608.7832, -53.96347, -208.59384],
        [6.93474, 4.6308937, 1.6512758],
        [8.253552, 5.426136, 0.9558364],
        [50.62288, 343.65768, -556.3136],
        [-539.3237, 343.65768, -115.611725],
        [7.073659, 5.426136, 1.8372401],
    ]
    .map(|arr| Vector::new(arr[0], arr[1], arr[2]));
    let convex_polyhedron_b =
        ConvexPolyhedron::from_convex_hull(&convex_polyhedron_b_points).unwrap();

    let num_contained_points = convex_polyhedron_a_points
        .iter()
        .filter(|point| {
            convex_polyhedron_b
                .project_local_point(**point, false)
                .is_inside
        })
        .count();

    let contact = query::contact(
        &Pose::IDENTITY,
        &convex_polyhedron_a,
        &Pose::IDENTITY,
        &convex_polyhedron_b,
        0.0,
    )
    .unwrap();

    assert!(num_contained_points == 1);
    assert!(contact.is_some());
}

// Shapes are pulled directly from the view frustums of a game engine.
#[test]
fn convex_polyhedra_intersection() {
    let convex_polyhedron_a_points = [
        [0.45838028, 5.7372417, 0.61019015],
        [1000.35846, -994.1628, 1000.51013],
        [1000.35834, -994.1628, -999.48987],
        [0.45838028, 5.7372417, 0.41019014],
        [0.45838028, 5.9372416, 0.61019015],
        [1000.35846, 1005.8372, 1000.51013],
        [1000.35834, 1005.8372, -999.48987],
        [0.45838028, 5.9372416, 0.41019014],
    ]
    .map(|arr| Vector::new(arr[0], arr[1], arr[2]));
    let convex_polyhedron_a =
        ConvexPolyhedron::from_convex_hull(&convex_polyhedron_a_points).unwrap();

    let convex_polyhedron_b_points = [
        [-0.4522132, 8.780189, 15.508082],
        [72537.12, -72390.1, -81439.88],
        [-74725.21, -72390.1, -79438.195],
        [-0.45368585, 8.780189, 15.508101],
        [-0.45221695, 8.780969, 15.507806],
        [72161.625, 5710.25, -109064.305],
        [-75100.7, 5710.25, -107062.62],
        [-0.4536896, 8.780969, 15.507825],
    ]
    .map(|arr| Vector::new(arr[0], arr[1], arr[2]));
    let convex_polyhedron_b =
        ConvexPolyhedron::from_convex_hull(&convex_polyhedron_b_points).unwrap();

    let num_contained_points = convex_polyhedron_a_points
        .iter()
        .filter(|point| {
            convex_polyhedron_b
                .project_local_point(**point, false)
                .is_inside
        })
        .count();

    let intersects = query::intersection_test(
        &Pose::IDENTITY,
        &convex_polyhedron_a,
        &Pose::IDENTITY,
        &convex_polyhedron_b,
    )
    .unwrap();

    assert!(num_contained_points == 4);
    assert!(intersects.intersecting);
}
