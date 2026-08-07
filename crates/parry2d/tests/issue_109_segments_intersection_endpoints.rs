// Regression test for https://github.com/dimforge/parry/issues/109
//
// `segments_intersection2d` used to classify exact segment endpoints as `OnEdge`
// instead of `OnVertex`, in both the collinear and the non-parallel branches.

use parry2d::math::Vector;
use parry2d::shape::SegmentPointLocation;
use parry2d::utils::{segments_intersection2d, SegmentsIntersection};

#[test]
fn identical_segments_intersect_on_vertices() {
    // Exact case from the issue.
    let a = Vector::new(10.0, 0.0);
    let b = Vector::new(10.0, 10.0);

    let intersection = segments_intersection2d(a, b, a, b, 0.0).unwrap();

    let SegmentsIntersection::Segment {
        first_loc1,
        first_loc2,
        second_loc1,
        second_loc2,
    } = intersection
    else {
        panic!("the intersection should be a Segment intersection");
    };

    assert_eq!(first_loc1, SegmentPointLocation::OnVertex(0));
    assert_eq!(first_loc2, SegmentPointLocation::OnVertex(0));
    assert_eq!(second_loc1, SegmentPointLocation::OnVertex(1));
    assert_eq!(second_loc2, SegmentPointLocation::OnVertex(1));
}

#[test]
fn identical_horizontal_segments_intersect_on_vertices() {
    // Same check through the non-vertical code path of `between()`.
    let a = Vector::new(-3.0, 2.0);
    let b = Vector::new(5.0, 2.0);

    let intersection = segments_intersection2d(a, b, a, b, 0.0).unwrap();

    let SegmentsIntersection::Segment {
        first_loc1,
        first_loc2,
        second_loc1,
        second_loc2,
    } = intersection
    else {
        panic!("the intersection should be a Segment intersection");
    };

    assert_eq!(first_loc1, SegmentPointLocation::OnVertex(0));
    assert_eq!(first_loc2, SegmentPointLocation::OnVertex(0));
    assert_eq!(second_loc1, SegmentPointLocation::OnVertex(1));
    assert_eq!(second_loc2, SegmentPointLocation::OnVertex(1));
}

#[test]
fn crossing_at_endpoint_is_on_vertex() {
    // Two non-parallel segments touching at seg1's vertex 1 == seg2's vertex 0.
    // Exercises the `s == 1.0` fix (previously `s == denom`, so `OnVertex(1)`
    // was unreachable).
    let a = Vector::new(0.0, 0.0);
    let b = Vector::new(1.0, 1.0);
    let c = Vector::new(1.0, 1.0);
    let d = Vector::new(2.0, 0.0);

    let intersection = segments_intersection2d(a, b, c, d, 0.0).unwrap();

    let SegmentsIntersection::Point { loc1, loc2 } = intersection else {
        panic!("the intersection should be a Point intersection");
    };

    assert_eq!(loc1, SegmentPointLocation::OnVertex(1));
    assert_eq!(loc2, SegmentPointLocation::OnVertex(0));
}

#[test]
fn crossing_at_interior_point_is_on_edge() {
    // Proper crossing in both interiors: still `OnEdge`.
    let a = Vector::new(-1.0, 0.0);
    let b = Vector::new(1.0, 0.0);
    let c = Vector::new(0.0, -1.0);
    let d = Vector::new(0.0, 1.0);

    let intersection = segments_intersection2d(a, b, c, d, 0.0).unwrap();

    let SegmentsIntersection::Point { loc1, loc2 } = intersection else {
        panic!("the intersection should be a Point intersection");
    };

    assert_eq!(loc1, SegmentPointLocation::OnEdge([0.5, 0.5]));
    assert_eq!(loc2, SegmentPointLocation::OnEdge([0.5, 0.5]));
}
