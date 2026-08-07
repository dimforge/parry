// Regression test for https://github.com/dimforge/parry/issues/76
//
// Point projection on a degenerate triangle used to skip every edge Voronoï test (the
// face normal being zero) and report the point as inside with distance 0; it now
// projects on the longest edge. Collinear cases adapted from PR #358.

use parry3d::math::Vector;
use parry3d::query::{PointQuery, PointQueryWithLocation};
use parry3d::shape::{Polyline, Shape, Triangle};

#[test]
fn degenerate_triangle_distance_matches_segment() {
    // Exact values from the issue.
    let p = Vector::new(1.10000002, -7.9000001, 16.5879993);

    let a = Vector::new(2.27699995, -7.9000001, 16.3180008);
    let b = Vector::new(-0.569999993, -8.10000038, 16.6070004);
    let c = Vector::new(-0.569999993, -8.10000038, 16.6070004);

    let line = Polyline::new(vec![a, b], None);
    let tri = Triangle::new(a, b, c);

    assert_eq!(tri.area(), 0.0);
    assert!(tri.compute_local_aabb().contains_local_point(p));

    let tri_dist = tri.distance_to_local_point(p, true);
    let line_dist = line.distance_to_local_point(p, true);

    assert!(line_dist > 0.0);
    assert_eq!(tri_dist, line_dist);

    // The non-solid projection must match too.
    assert_eq!(
        tri.distance_to_local_point(p, false),
        line.distance_to_local_point(p, false)
    );
}

#[test]
fn degenerate_triangle_distance_several_points() {
    let a = Vector::new(-1.0, 2.0, 0.5);
    let b = Vector::new(3.0, -1.0, 2.0);

    // b == c: the degenerate triangle must behave like its longest segment.
    let tri = Triangle::new(a, b, b);
    let seg = parry3d::shape::Segment::new(a, b);

    let queries = [
        Vector::new(0.0, 0.0, 0.0),
        Vector::new(10.0, 10.0, -3.0),
        Vector::new(-5.0, 2.5, 1.0),
        Vector::new(1.0, 0.5, 1.25), // Near the middle of the segment.
        a,                           // Exactly on a vertex.
    ];

    for pt in queries {
        let tri_dist = tri.distance_to_local_point(pt, true);
        let seg_dist = seg.distance_to_local_point(pt, true);
        assert!(
            (tri_dist - seg_dist).abs() <= 1.0e-6,
            "point {pt:?}: triangle dist {tri_dist} != segment dist {seg_dist}"
        );
        assert!(tri_dist.is_finite());
    }
}

// Cases from #358 by its author.
#[test]
fn two_identical_points_triangle_projection_is_finite() {
    let triangle = Triangle::new(
        Vector::new(40.0, 0.0, 0.0),
        Vector::new(0.0, 80.0, 0.0),
        Vector::new(0.0, 80.0, 0.0),
    );

    let res = triangle.project_local_point_and_get_location(Vector::new(10.0, 20.0, 0.0), false);
    assert!(res.0.point.is_finite());

    let res = triangle.project_local_point_and_get_location(Vector::new(40.0, 0.0, 0.0), false);
    assert!(res.0.point.is_finite());
}

// Cases from #358 by its author.
#[test]
fn collinear_points_triangle_projection_is_finite() {
    let triangle = Triangle::new(
        Vector::new(0.0, 0.0, 0.0),
        Vector::new(100.0, 0.0, 0.0),
        Vector::new(160.0, 0.0, 0.0),
    );

    // Point considered "inside" (on the line).
    let res = triangle.project_local_point_and_get_location(Vector::new(10.0, 0.0, 0.0), false);
    assert!(res.0.is_inside);
    assert!(res.0.point.is_finite());

    // Point off the line: not inside, finite projection.
    let res = triangle.project_local_point_and_get_location(Vector::new(10.0, 10.0, 0.0), false);
    assert!(!res.0.is_inside);
    assert!(res.0.point.is_finite());

    // The solid flag must not change anything for a degenerate triangle:
    // it has no interior.
    let res = triangle.project_local_point_and_get_location(Vector::new(10.0, 10.0, 0.0), true);
    assert!(!res.0.is_inside);
    assert_eq!(
        triangle.distance_to_local_point(Vector::new(10.0, 10.0, 0.0), true),
        10.0
    );
}
