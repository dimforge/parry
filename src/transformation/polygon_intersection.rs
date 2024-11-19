use log::error;
use na::Point2;
use ordered_float::OrderedFloat;

use crate::math::Real;
use crate::shape::{SegmentPointLocation, Triangle, TriangleOrientation};
use crate::utils::hashmap::HashMap;
use crate::utils::{self, SegmentsIntersection};

#[derive(Copy, Clone, PartialEq, Debug)]
pub struct PolygonIntersectionTolerances {
    /// The epsilon given to [`Triangle::orientation2d`] for detecting if three points are collinear.
    ///
    /// Three points forming a triangle with an area smaller than this value are considered collinear.
    pub collinearity_epsilon: Real,
}

impl Default for PolygonIntersectionTolerances {
    fn default() -> Self {
        Self {
            collinearity_epsilon: Real::EPSILON * 100.0,
        }
    }
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
enum InFlag {
    // The current neighborhood of the traversed point on poly1 is inside poly2.
    Poly1IsInside,
    // The current neighborhood of the traversed point on poly2 is inside poly1.
    Poly2IsInside,
    Unknown,
}

/// Location of a point on a polyline.
#[derive(Copy, Clone, Debug, PartialEq)]
pub enum PolylinePointLocation {
    /// Point on a vertex.
    OnVertex(usize),
    /// Point on an edge.
    OnEdge(usize, usize, [Real; 2]),
}

impl PolylinePointLocation {
    /// The barycentric coordinates such that the point in the intersected segment `[a, b]` is
    /// equal to `a + (b - a) * centered_bcoords`.
    fn centered_bcoords(&self, edge: [usize; 2]) -> Real {
        match self {
            Self::OnVertex(vid) => {
                if *vid == edge[0] {
                    0.0
                } else {
                    1.0
                }
            }
            Self::OnEdge(ia, ib, bcoords) => {
                assert_eq!([*ia, *ib], edge);
                bcoords[1]
            }
        }
    }

    /// Computes the point corresponding to this location.
    pub fn to_point(self, pts: &[Point2<Real>]) -> Point2<Real> {
        match self {
            PolylinePointLocation::OnVertex(i) => pts[i],
            PolylinePointLocation::OnEdge(i1, i2, bcoords) => {
                pts[i1] * bcoords[0] + pts[i2].coords * bcoords[1]
            }
        }
    }

    fn from_segment_point_location(a: usize, b: usize, loc: SegmentPointLocation) -> Self {
        match loc {
            SegmentPointLocation::OnVertex(0) => PolylinePointLocation::OnVertex(a),
            SegmentPointLocation::OnVertex(1) => PolylinePointLocation::OnVertex(b),
            SegmentPointLocation::OnVertex(_) => unreachable!(),
            SegmentPointLocation::OnEdge(bcoords) => PolylinePointLocation::OnEdge(a, b, bcoords),
        }
    }
}

/// Computes the intersection points of two convex polygons.
///
/// The resulting polygon is output vertex-by-vertex to the `out` vector.
/// This is the same as [`convex_polygons_intersection_points_with_tolerances`] with the tolerances
/// set to their default values.
pub fn convex_polygons_intersection_points(
    poly1: &[Point2<Real>],
    poly2: &[Point2<Real>],
    out: &mut Vec<Point2<Real>>,
) {
    convex_polygons_intersection_points_with_tolerances(poly1, poly2, Default::default(), out);
}

/// Computes the intersection points of two convex polygons.
///
/// The resulting polygon is output vertex-by-vertex to the `out` vector.
pub fn convex_polygons_intersection_points_with_tolerances(
    poly1: &[Point2<Real>],
    poly2: &[Point2<Real>],
    tolerances: PolygonIntersectionTolerances,
    out: &mut Vec<Point2<Real>>,
) {
    convex_polygons_intersection_with_tolerances(poly1, poly2, tolerances, |loc1, loc2| {
        if let Some(loc1) = loc1 {
            out.push(loc1.to_point(poly1))
        } else if let Some(loc2) = loc2 {
            out.push(loc2.to_point(poly2))
        }
    })
}

/// Computes the intersection of two convex polygons.
///
/// The resulting polygon is output vertex-by-vertex to the `out` closure.
/// This is the same as [`convex_polygons_intersection_with_tolerances`] with the tolerances
/// set to their default values.
pub fn convex_polygons_intersection(
    poly1: &[Point2<Real>],
    poly2: &[Point2<Real>],
    out: impl FnMut(Option<PolylinePointLocation>, Option<PolylinePointLocation>),
) {
    convex_polygons_intersection_with_tolerances(poly1, poly2, Default::default(), out)
}

/// Computes the intersection of two convex polygons.
///
/// The resulting polygon is output vertex-by-vertex to the `out` closure.
pub fn convex_polygons_intersection_with_tolerances(
    poly1: &[Point2<Real>],
    poly2: &[Point2<Real>],
    tolerances: PolygonIntersectionTolerances,
    mut out: impl FnMut(Option<PolylinePointLocation>, Option<PolylinePointLocation>),
) {
    // TODO: this does not handle correctly the case where the
    // first triangle of the polygon is degenerate.
    let rev1 = poly1.len() > 2
        && Triangle::orientation2d(
            &poly1[0],
            &poly1[1],
            &poly1[2],
            tolerances.collinearity_epsilon,
        ) == TriangleOrientation::Clockwise;
    let rev2 = poly2.len() > 2
        && Triangle::orientation2d(
            &poly2[0],
            &poly2[1],
            &poly2[2],
            tolerances.collinearity_epsilon,
        ) == TriangleOrientation::Clockwise;

    let len1 = poly1.len();
    let len2 = poly2.len();

    let mut i1 = 0; // Current index on the first polyline.
    let mut i2 = 0; // Current index on the second polyline.
    let mut nsteps1 = 0; // Number of times we advanced on the first polyline.
    let mut nsteps2 = 0; // Number of times we advanced on the second polyline.
    let mut inflag = InFlag::Unknown;
    let mut first_point_found = false;

    // Quit when both adv. indices have cycled, or one has cycled twice.
    while (nsteps1 < len1 || nsteps2 < len2) && nsteps1 < 2 * len1 && nsteps2 < 2 * len2 {
        let (a1, b1) = if rev1 {
            ((len1 - i1) % len1, len1 - i1 - 1)
        } else {
            // Point before `i1`, and point at `i1`.
            ((i1 + len1 - 1) % len1, i1)
        };

        let (a2, b2) = if rev2 {
            ((len2 - i2) % len2, len2 - i2 - 1)
        } else {
            // Point before `i2`, and point at `i2`.
            ((i2 + len2 - 1) % len2, i2)
        };

        let dir_edge1 = poly1[b1] - poly1[a1];
        let dir_edge2 = poly2[b2] - poly2[a2];

        // If there is an intersection, this will determine if the edge from poly2 is transitioning
        // Left -> Right (CounterClockwise) or Right -> Left (Clockwise) relative to the edge from
        // poly1.
        let cross = Triangle::orientation2d(
            &Point2::origin(),
            &Point2::from(dir_edge1),
            &Point2::from(dir_edge2),
            tolerances.collinearity_epsilon,
        );
        // Determines if b1 is left (CounterClockwise) or right (Clockwise) of [a2, b2].
        let a2_b2_b1 = Triangle::orientation2d(
            &poly2[a2],
            &poly2[b2],
            &poly1[b1],
            tolerances.collinearity_epsilon,
        );
        // Determines if b2 is left (CounterClockwise) or right (Clockwise) of [a1, b1].
        let a1_b1_b2 = Triangle::orientation2d(
            &poly1[a1],
            &poly1[b1],
            &poly2[b2],
            tolerances.collinearity_epsilon,
        );

        // If edge1 & edge2 intersect, update inflag.
        if let Some(inter) = utils::segments_intersection2d(
            &poly1[a1],
            &poly1[b1],
            &poly2[a2],
            &poly2[b2],
            tolerances.collinearity_epsilon,
        ) {
            match inter {
                SegmentsIntersection::Point { loc1, loc2 } => {
                    if a2_b2_b1 != TriangleOrientation::Degenerate
                        && a1_b1_b2 != TriangleOrientation::Degenerate
                    {
                        let loc1 = PolylinePointLocation::from_segment_point_location(a1, b1, loc1);
                        let loc2 = PolylinePointLocation::from_segment_point_location(a2, b2, loc2);
                        out(Some(loc1), Some(loc2));

                        if inflag == InFlag::Unknown && !first_point_found {
                            // This is the first point, reset the number of steps since we are
                            // effectively starting the actual traversal now.
                            nsteps1 = 0;
                            nsteps2 = 0;
                            first_point_found = true;
                        }

                        if a2_b2_b1 == TriangleOrientation::CounterClockwise {
                            // The point b1 is left of [a2, b2] so it is inside poly2 ???
                            inflag = InFlag::Poly1IsInside;
                        } else if a1_b1_b2 == TriangleOrientation::CounterClockwise {
                            // The point b2 is left of [a1, b1] so it is inside poly1 ???
                            inflag = InFlag::Poly2IsInside;
                        }
                    }
                }
                SegmentsIntersection::Segment {
                    first_loc1,
                    first_loc2,
                    second_loc1,
                    second_loc2,
                } => {
                    if dir_edge1.dot(&dir_edge2) < 0.0 {
                        // Special case: edge1 & edge2 overlap and oppositely oriented. The
                        //               intersection is degenerate (equals to a segment). Output
                        //               the segment and exit.
                        let loc1 =
                            PolylinePointLocation::from_segment_point_location(a1, b1, first_loc1);
                        let loc2 =
                            PolylinePointLocation::from_segment_point_location(a2, b2, first_loc2);
                        out(Some(loc1), Some(loc2));

                        let loc1 =
                            PolylinePointLocation::from_segment_point_location(a1, b1, second_loc1);
                        let loc2 =
                            PolylinePointLocation::from_segment_point_location(a2, b2, second_loc2);
                        out(Some(loc1), Some(loc2));
                        return;
                    }
                }
            }
        }

        // Special case: edge1 & edge2 parallel and separated.
        if cross == TriangleOrientation::Degenerate
            && a2_b2_b1 == TriangleOrientation::Clockwise
            && a1_b1_b2 == TriangleOrientation::Clockwise
        // TODO: should this also include any case where a2_b2_b1 and a1_b1_b2 are both different from Degenerate?
        {
            return;
        }
        // Special case: edge1 & edge2 collinear.
        else if cross == TriangleOrientation::Degenerate
            && a2_b2_b1 == TriangleOrientation::Degenerate
            && a1_b1_b2 == TriangleOrientation::Degenerate
        {
            // Advance but do not output point.
            if inflag == InFlag::Poly1IsInside {
                i2 = advance(i2, &mut nsteps2, len2);
            } else {
                i1 = advance(i1, &mut nsteps1, len1);
            }
        }
        // Generic cases.
        else if cross == TriangleOrientation::CounterClockwise {
            if a1_b1_b2 == TriangleOrientation::CounterClockwise {
                if inflag == InFlag::Poly1IsInside {
                    out(Some(PolylinePointLocation::OnVertex(b1)), None)
                }
                i1 = advance(i1, &mut nsteps1, len1);
            } else {
                if inflag == InFlag::Poly2IsInside {
                    out(None, Some(PolylinePointLocation::OnVertex(b2)))
                }
                i2 = advance(i2, &mut nsteps2, len2);
            }
        } else {
            // We have cross == TriangleOrientation::Clockwise.
            if a2_b2_b1 == TriangleOrientation::CounterClockwise {
                if inflag == InFlag::Poly2IsInside {
                    out(None, Some(PolylinePointLocation::OnVertex(b2)))
                }
                i2 = advance(i2, &mut nsteps2, len2);
            } else {
                if inflag == InFlag::Poly1IsInside {
                    out(Some(PolylinePointLocation::OnVertex(b1)), None)
                }
                i1 = advance(i1, &mut nsteps1, len1);
            }
        }
    }

    if !first_point_found {
        // No intersection: test if one polygon completely encloses the other.
        let mut orient = TriangleOrientation::Degenerate;
        let mut ok = true;

        // TODO: avoid the O(n²) complexity.
        for a in 0..len1 {
            for p2 in poly2 {
                let a_minus_1 = (a + len1 - 1) % len1; // a - 1
                let new_orient = Triangle::orientation2d(
                    &poly1[a_minus_1],
                    &poly1[a],
                    p2,
                    tolerances.collinearity_epsilon,
                );

                if orient == TriangleOrientation::Degenerate {
                    orient = new_orient
                } else if new_orient != orient && new_orient != TriangleOrientation::Degenerate {
                    ok = false;
                    break;
                }
            }
        }

        if ok {
            for mut b in 0..len2 {
                if rev2 {
                    b = len2 - b - 1;
                }
                out(None, Some(PolylinePointLocation::OnVertex(b)))
            }
        }

        let mut orient = TriangleOrientation::Degenerate;
        let mut ok = true;

        // TODO: avoid the O(n²) complexity.
        for b in 0..len2 {
            for p1 in poly1 {
                let b_minus_1 = (b + len2 - 1) % len2; // = b - 1
                let new_orient = Triangle::orientation2d(
                    &poly2[b_minus_1],
                    &poly2[b],
                    p1,
                    tolerances.collinearity_epsilon,
                );

                if orient == TriangleOrientation::Degenerate {
                    orient = new_orient
                } else if new_orient != orient && new_orient != TriangleOrientation::Degenerate {
                    ok = false;
                    break;
                }
            }
        }

        if ok {
            for mut a in 0..len1 {
                if rev1 {
                    a = len1 - a - 1;
                }
                out(Some(PolylinePointLocation::OnVertex(a)), None)
            }
        }
    }
}

#[inline]
fn advance(a: usize, aa: &mut usize, n: usize) -> usize {
    *aa += 1;
    (a + 1) % n
}

#[derive(thiserror::Error, Debug)]
pub enum PolygonsIntersectionError {
    #[error("Infinite loop detected; input polygons are ill-formed.")]
    InfiniteLoop,
}

/// Compute intersections between two polygons that may be non-convex but that must not self-intersect.
///
/// The input polygons are assumed to not self-intersect, and to be oriented counter-clockwise.
///
/// The resulting polygon is output vertex-by-vertex to the `out` closure.
/// If two `None` are given to the `out` closure, then one connected component of the intersection
/// polygon is complete.
///
/// If the polygons are known to be convex, use [`convex_polygons_intersection_points`] instead for better
/// performances.
pub fn polygons_intersection_points(
    poly1: &[Point2<Real>],
    poly2: &[Point2<Real>],
) -> Result<Vec<Vec<Point2<Real>>>, PolygonsIntersectionError> {
    let mut result = vec![];
    let mut curr_poly = vec![];
    polygons_intersection(poly1, poly2, |loc1, loc2| {
        if let Some(loc1) = loc1 {
            curr_poly.push(loc1.to_point(poly1))
        } else if let Some(loc2) = loc2 {
            curr_poly.push(loc2.to_point(poly2))
        } else if !curr_poly.is_empty() {
            result.push(std::mem::take(&mut curr_poly));
        }
    })?;

    Ok(result)
}

/// Compute intersections between two polygons that may be non-convex but that must not self-intersect.
///
/// The input polygons are assumed to not self-intersect, and to be oriented counter-clockwise.
///
/// The resulting polygon is output vertex-by-vertex to the `out` closure.
/// If two `None` are given to the `out` closure, then one connected component of the intersection
/// polygon is complete.
///
/// If the polygons are known to be convex, use [`convex_polygons_intersection`] instead for better
/// performances.
pub fn polygons_intersection(
    poly1: &[Point2<Real>],
    poly2: &[Point2<Real>],
    mut out: impl FnMut(Option<PolylinePointLocation>, Option<PolylinePointLocation>),
) -> Result<(), PolygonsIntersectionError> {
    let tolerances = PolygonIntersectionTolerances::default();

    #[derive(Debug)]
    struct ToTraverse {
        poly: usize,
        edge: EdgeId,
    }

    let (intersections, num_intersections) =
        compute_sorted_edge_intersections(poly1, poly2, tolerances.collinearity_epsilon);
    let mut visited = vec![false; num_intersections];
    let segment = |eid: EdgeId, poly: &[Point2<Real>]| [poly[eid], poly[(eid + 1) % poly.len()]];

    // Traverse all the intersections.
    for inters in intersections[0].values() {
        for inter in inters {
            if visited[inter.id] {
                continue;
            }

            // We found an intersection we haven’t visited yet, traverse the loop, alternating
            // between poly1 and poly2 when reaching an intersection.
            let [a1, b1] = segment(inter.edges[0], poly1);
            let [a2, b2] = segment(inter.edges[1], poly2);
            let poly_to_traverse =
                match Triangle::orientation2d(&a1, &b1, &a2, tolerances.collinearity_epsilon) {
                    TriangleOrientation::Clockwise => 1,
                    TriangleOrientation::CounterClockwise => 0,
                    TriangleOrientation::Degenerate => {
                        match Triangle::orientation2d(
                            &a1,
                            &b1,
                            &b2,
                            tolerances.collinearity_epsilon,
                        ) {
                            TriangleOrientation::Clockwise => 0,
                            TriangleOrientation::CounterClockwise => 1,
                            TriangleOrientation::Degenerate => {
                                log::debug!("Unhandled edge-edge overlap case.");
                                0
                            }
                        }
                    }
                };

            #[derive(Debug)]
            enum TraversalStatus {
                OnVertex,
                OnIntersection(usize),
            }

            let polys = [poly1, poly2];
            let mut to_traverse = ToTraverse {
                poly: poly_to_traverse,
                edge: inter.edges[poly_to_traverse],
            };

            let mut status = TraversalStatus::OnIntersection(inter.id);

            for loop_id in 0.. {
                if loop_id > poly1.len() * poly2.len() {
                    return Err(PolygonsIntersectionError::InfiniteLoop);
                }

                let empty = Vec::new();
                let edge_inters = intersections[to_traverse.poly]
                    .get(&to_traverse.edge)
                    .unwrap_or(&empty);

                match status {
                    TraversalStatus::OnIntersection(inter_id) => {
                        let (curr_inter_pos, curr_inter) = edge_inters
                            .iter()
                            .enumerate()
                            .find(|(_, inter)| inter.id == inter_id)
                            .unwrap_or_else(|| unreachable!());

                        if visited[curr_inter.id] {
                            // We already saw this intersection: we looped back to the start of
                            // the intersection polygon.
                            out(None, None);
                            break;
                        }

                        out(Some(curr_inter.locs[0]), Some(curr_inter.locs[1]));
                        visited[curr_inter.id] = true;

                        if curr_inter_pos + 1 < edge_inters.len() {
                            // There are other intersections after this one.
                            // Move forward to the next intersection point and move on to traversing
                            // the other polygon.
                            let next_inter = &edge_inters[curr_inter_pos + 1];
                            to_traverse.poly = (to_traverse.poly + 1) % 2;
                            to_traverse.edge = next_inter.edges[to_traverse.poly];
                            status = TraversalStatus::OnIntersection(next_inter.id);
                        } else {
                            // This was the last intersection, move to the next vertex on the
                            // same polygon.
                            to_traverse.edge =
                                (to_traverse.edge + 1) % polys[to_traverse.poly].len();
                            status = TraversalStatus::OnVertex;
                        }
                    }
                    TraversalStatus::OnVertex => {
                        let location = PolylinePointLocation::OnVertex(to_traverse.edge);

                        if to_traverse.poly == 0 {
                            out(Some(location), None);
                        } else {
                            out(None, Some(location))
                        };

                        if let Some(first_intersection) = edge_inters.first() {
                            // Jump on the first intersection and move on to the other polygon.
                            to_traverse.poly = (to_traverse.poly + 1) % 2;
                            to_traverse.edge = first_intersection.edges[to_traverse.poly];
                            status = TraversalStatus::OnIntersection(first_intersection.id);
                        } else {
                            // Move forward to the next vertex/edge on the same polygon.
                            to_traverse.edge =
                                (to_traverse.edge + 1) % polys[to_traverse.poly].len();
                        }
                    }
                }
            }
        }
    }

    // If there are no intersection, check if one polygon is inside the other.
    if intersections[0].is_empty() {
        if utils::point_in_poly2d(&poly1[0], poly2) {
            for pt_id in 0..poly1.len() {
                out(Some(PolylinePointLocation::OnVertex(pt_id)), None)
            }
            out(None, None);
        } else if utils::point_in_poly2d(&poly2[0], poly1) {
            for pt_id in 0..poly2.len() {
                out(None, Some(PolylinePointLocation::OnVertex(pt_id)))
            }
            out(None, None);
        }
    }

    Ok(())
}

type EdgeId = usize;
type IntersectionId = usize;

#[derive(Copy, Clone, Debug)]
struct IntersectionPoint {
    id: IntersectionId,
    edges: [EdgeId; 2],
    locs: [PolylinePointLocation; 2],
}

fn compute_sorted_edge_intersections(
    poly1: &[Point2<Real>],
    poly2: &[Point2<Real>],
    eps: Real,
) -> ([HashMap<EdgeId, Vec<IntersectionPoint>>; 2], usize) {
    let mut inter1: HashMap<EdgeId, Vec<IntersectionPoint>> = HashMap::default();
    let mut inter2: HashMap<EdgeId, Vec<IntersectionPoint>> = HashMap::default();
    let mut id = 0;

    // Find the intersections.
    // TODO: this is a naive O(n²) check. Could use an acceleration structure for large polygons.
    for i1 in 0..poly1.len() {
        let j1 = (i1 + 1) % poly1.len();

        for i2 in 0..poly2.len() {
            let j2 = (i2 + 1) % poly2.len();

            let Some(inter) =
                utils::segments_intersection2d(&poly1[i1], &poly1[j1], &poly2[i2], &poly2[j2], eps)
            else {
                continue;
            };

            match inter {
                SegmentsIntersection::Point { loc1, loc2 } => {
                    let loc1 = PolylinePointLocation::from_segment_point_location(i1, j1, loc1);
                    let loc2 = PolylinePointLocation::from_segment_point_location(i2, j2, loc2);
                    let intersection = IntersectionPoint {
                        id,
                        edges: [i1, i2],
                        locs: [loc1, loc2],
                    };
                    inter1.entry(i1).or_default().push(intersection);
                    inter2.entry(i2).or_default().push(intersection);
                    id += 1;
                }
                SegmentsIntersection::Segment { .. } => {
                    // TODO
                    log::debug!(
                        "Collinear segment-segment intersections not properly handled yet."
                    );
                }
            }
        }
    }

    // Sort the intersections.
    for inters in inter1.values_mut() {
        inters.sort_by_key(|a| {
            let edge = [a.edges[0], (a.edges[0] + 1) % poly1.len()];
            OrderedFloat(a.locs[0].centered_bcoords(edge))
        });
    }

    for inters in inter2.values_mut() {
        inters.sort_by_key(|a| {
            let edge = [a.edges[1], (a.edges[1] + 1) % poly2.len()];
            OrderedFloat(a.locs[1].centered_bcoords(edge))
        });
    }

    ([inter1, inter2], id)
}

#[cfg(all(test, feature = "dim2"))]
mod test {
    use crate::query::PointQuery;
    use crate::shape::Triangle;
    use crate::transformation::convex_polygons_intersection_points_with_tolerances;
    use crate::transformation::polygon_intersection::PolygonIntersectionTolerances;
    use na::Point2;

    #[test]
    fn intersect_triangle_common_vertex() {
        let tris = [
            (
                Triangle::new(
                    Point2::new(-0.0008759537858568062, -2.0103871966663305),
                    Point2::new(0.3903908709629763, -1.3421764825890266),
                    Point2::new(1.3380817875388151, -2.0098007857739013),
                ),
                Triangle::new(
                    Point2::new(0.0, -0.0),
                    Point2::new(-0.0008759537858568062, -2.0103871966663305),
                    Point2::new(1.9991979155226394, -2.009511242880474),
                ),
            ),
            (
                Triangle::new(
                    Point2::new(0.7319315811016305, -0.00004046981523721891),
                    Point2::new(2.0004914907008944, -0.00011061077714557787),
                    Point2::new(1.1848406021956144, -0.8155712451545468),
                ),
                Triangle::new(
                    Point2::new(0.0, 0.0),
                    Point2::new(0.00011061077714557787, -2.000024893134292),
                    Point2::new(2.0004914907008944, -0.00011061077714557787),
                ),
            ),
            (
                Triangle::new(
                    Point2::new(-0.000049995168258705205, -0.9898801451981707),
                    Point2::new(0.0, -0.0),
                    Point2::new(0.583013294019752, -1.4170136900568633),
                ),
                Triangle::new(
                    Point2::new(0.0, -0.0),
                    Point2::new(-0.00010101395240669591, -2.000027389553894),
                    Point2::new(2.000372544168497, 0.00010101395240669591),
                ),
            ),
            (
                Triangle::new(
                    Point2::new(-0.940565646581769, -0.939804943675256),
                    Point2::new(0.0, -0.0),
                    Point2::new(-0.001533592366792066, -0.9283586484736431),
                ),
                Triangle::new(
                    Point2::new(0.0, -0.0),
                    Point2::new(-2.00752629829582, -2.0059026672784825),
                    Point2::new(-0.0033081650580626698, -2.0025945022204197),
                ),
            ),
        ];

        for (tri1, tri2) in tris {
            let mut inter = Vec::new();
            let tolerances = PolygonIntersectionTolerances {
                collinearity_epsilon: 1.0e-5,
            };
            convex_polygons_intersection_points_with_tolerances(
                tri1.vertices(),
                tri2.vertices(),
                tolerances,
                &mut inter,
            );

            println!("----------");
            println!("inter: {:?}", inter);
            println!(
                "tri1 is in tri2: {}",
                tri1.vertices().iter().all(|pt| tri2
                    .project_local_point(pt, false)
                    .is_inside_eps(pt, 1.0e-5))
            );
            println!(
                "tri2 is in tri1: {}",
                tri2.vertices().iter().all(|pt| tri1
                    .project_local_point(pt, false)
                    .is_inside_eps(pt, 1.0e-5))
            );
            for pt in &inter {
                let proj1 = tri1.project_local_point(&pt, false);
                let proj2 = tri2.project_local_point(&pt, false);
                assert!(proj1.is_inside_eps(&pt, 1.0e-5));
                assert!(proj2.is_inside_eps(&pt, 1.0e-5));
            }
        }
    }
}
