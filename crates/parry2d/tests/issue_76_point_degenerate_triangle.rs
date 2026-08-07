// Regression test for https://github.com/dimforge/parry/issues/76 (2D version).
//
// A degenerate triangle reaching the face Voronoï region used to be handled by the
// `solid` branch, reporting the point as inside with distance 0 (or producing NaNs);
// it now projects on the longest edge. Collinear cases adapted from PR #358.

use parry2d::math::Vector;
use parry2d::query::{PointQuery, PointQueryWithLocation};
use parry2d::shape::{Segment, Triangle};

#[test]
fn degenerate_triangle_distance_matches_segment() {
    let a = Vector::new(-1.0, 2.0);
    let b = Vector::new(3.0, -1.0);

    // b == c: the degenerate triangle must behave like its longest segment.
    let tri = Triangle::new(a, b, b);
    let seg = Segment::new(a, b);

    let queries = [
        Vector::new(0.0, 0.0),
        Vector::new(10.0, 10.0),
        Vector::new(-5.0, 2.5),
        Vector::new(1.0, 0.5), // Near the middle of the segment.
        a,                     // Exactly on a vertex.
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
        Vector::new(40.0, 0.0),
        Vector::new(0.0, 80.0),
        Vector::new(0.0, 80.0),
    );

    let res = triangle.project_local_point_and_get_location(Vector::new(10.0, 20.0), false);
    assert!(res.0.point.is_finite());

    let res = triangle.project_local_point_and_get_location(Vector::new(40.0, 0.0), false);
    assert!(res.0.point.is_finite());
}

// Cases from #358 by its author.
#[test]
fn collinear_points_triangle_projection_is_finite() {
    let triangle = Triangle::new(
        Vector::new(0.0, 0.0),
        Vector::new(100.0, 0.0),
        Vector::new(160.0, 0.0),
    );

    // Point considered "inside" (on the line).
    let res = triangle.project_local_point_and_get_location(Vector::new(10.0, 0.0), false);
    assert!(res.0.is_inside);
    assert!(res.0.point.is_finite());

    // Point off the line: not inside, finite projection, correct distance
    // even with `solid = true` (a degenerate triangle has no interior).
    let res = triangle.project_local_point_and_get_location(Vector::new(10.0, 10.0), false);
    assert!(!res.0.is_inside);
    assert!(res.0.point.is_finite());
    assert_eq!(
        triangle.distance_to_local_point(Vector::new(10.0, 10.0), true),
        10.0
    );
}
