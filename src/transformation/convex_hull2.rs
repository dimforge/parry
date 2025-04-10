use alloc::vec::Vec;
use core::marker::PhantomData;

use crate::math::Real;
use crate::transformation::convex_hull_utils::{indexed_support_point_id, support_point_id};
use na::{self, Point2, Vector2};
use num_traits::Zero;

/// Computes the convex hull of a set of 2d points.
///
/// The computed convex-hull have its points given in counter-clockwise order.
#[cfg(feature = "dim2")]
pub fn convex_hull2(points: &[Point2<Real>]) -> Vec<Point2<Real>> {
    convex_hull2_idx(points)
        .into_iter()
        .map(|id| points[id])
        .collect()
}

/// Computes the convex hull of a set of 2d points and returns only the indices of the hull
/// vertices.
///
/// The computed convex-hull have its points given in counter-clockwise order.
pub fn convex_hull2_idx(points: &[Point2<Real>]) -> Vec<usize> {
    let mut undecidable_points = Vec::new();
    let mut segments = get_initial_polyline(points, &mut undecidable_points);

    let mut i = 0;
    while i != segments.len() {
        if !segments[i].valid {
            i += 1;
            continue;
        }

        let pt_id = indexed_support_point_id(
            &segments[i].normal,
            points,
            segments[i].visible_points.iter().copied(),
        );

        if let Some(point) = pt_id {
            segments[i].valid = false;

            attach_and_push_facets2(
                segments[i].prev,
                segments[i].next,
                point,
                points,
                &mut segments,
                i,
                &mut undecidable_points,
            );
        }

        i += 1;
    }

    let mut idx = Vec::new();
    let mut curr_facet = 0;

    while !segments[curr_facet].valid {
        curr_facet += 1
    }

    let first_facet = curr_facet;

    loop {
        let curr = &segments[curr_facet];

        if curr.valid {
            idx.push(curr.pts[0]);
        }

        curr_facet = curr.next;

        if curr_facet == first_facet {
            break;
        }
    }

    idx
}

fn get_initial_polyline(
    points: &[Point2<Real>],
    undecidable: &mut Vec<usize>,
) -> Vec<SegmentFacet> {
    let mut res = Vec::new();

    assert!(points.len() >= 2);

    let p1 = support_point_id(&Vector2::x(), points).unwrap();
    let mut p2 = p1;

    let direction = [-Vector2::x(), -Vector2::y(), Vector2::y()];

    for dir in direction.iter() {
        p2 = support_point_id(dir, points).unwrap();

        let p1p2 = points[p2] - points[p1];

        if !p1p2.norm_squared().is_zero() {
            break;
        }
    }

    assert!(
        p1 != p2,
        "Failed to build the 2d convex hull of this point cloud."
    );

    // Build two facets with opposite normals.
    let mut f1 = SegmentFacet::new(p1, p2, 1, 1, points);
    let mut f2 = SegmentFacet::new(p2, p1, 0, 0, points);

    // Attribute points to each facet.
    for i in 0..points.len() {
        if i == p1 || i == p2 {
            continue;
        }
        if f1.can_be_seen_by(i, points) {
            f1.visible_points.push(i);
        } else if f2.can_be_seen_by(i, points) {
            f2.visible_points.push(i);
        } else {
            // The point is collinear.
            undecidable.push(i);
        }
    }

    res.push(f1);
    res.push(f2);

    res
}

fn attach_and_push_facets2(
    prev_facet: usize,
    next_facet: usize,
    point: usize,
    points: &[Point2<Real>],
    segments: &mut Vec<SegmentFacet>,
    removed_facet: usize,
    undecidable: &mut Vec<usize>,
) {
    let new_facet1_id = segments.len();
    let new_facet2_id = new_facet1_id + 1;
    let prev_pt = segments[prev_facet].pts[1];
    let next_pt = segments[next_facet].pts[0];

    let mut new_facet1 = SegmentFacet::new(prev_pt, point, prev_facet, new_facet2_id, points);
    let mut new_facet2 = SegmentFacet::new(point, next_pt, new_facet1_id, next_facet, points);

    segments[prev_facet].next = new_facet1_id;
    segments[next_facet].prev = new_facet2_id;

    // Assign to each facets some of the points which can see it.
    for visible_point in segments[removed_facet].visible_points.iter() {
        if *visible_point == point {
            continue;
        }

        if new_facet1.can_be_seen_by(*visible_point, points) {
            new_facet1.visible_points.push(*visible_point);
        } else if new_facet2.can_be_seen_by(*visible_point, points) {
            new_facet2.visible_points.push(*visible_point);
        }
        // If none of the facet can be seen from the point, it is naturally deleted.
    }

    // Try to assign collinear points to one of the new facets
    let mut i = 0;

    while i != undecidable.len() {
        if new_facet1.can_be_seen_by(undecidable[i], points) {
            new_facet1.visible_points.push(undecidable[i]);
            let _ = undecidable.swap_remove(i);
        } else if new_facet2.can_be_seen_by(undecidable[i], points) {
            new_facet2.visible_points.push(undecidable[i]);
            let _ = undecidable.swap_remove(i);
        } else {
            i += 1;
        }
    }

    segments.push(new_facet1);
    segments.push(new_facet2);
}

struct SegmentFacet {
    pub valid: bool,
    pub normal: Vector2<Real>,
    pub next: usize,
    pub prev: usize,
    pub pts: [usize; 2],
    pub visible_points: Vec<usize>,
    pt_type: PhantomData<Point2<Real>>,
}

impl SegmentFacet {
    pub fn new(
        p1: usize,
        p2: usize,
        prev: usize,
        next: usize,
        points: &[Point2<Real>],
    ) -> SegmentFacet {
        let p1p2 = points[p2] - points[p1];

        let mut normal = Vector2::new(p1p2.y, -p1p2.x);
        let norm = normal.normalize_mut();

        SegmentFacet {
            valid: norm != 0.0,
            normal,
            prev,
            next,
            pts: [p1, p2],
            visible_points: Vec::new(),
            pt_type: PhantomData,
        }
    }

    pub fn can_be_seen_by(&self, point: usize, points: &[Point2<Real>]) -> bool {
        let p0 = &points[self.pts[0]];
        let pt = &points[point];

        let _eps = crate::math::DEFAULT_EPSILON;

        (*pt - *p0).dot(&self.normal) > _eps * na::convert::<f64, Real>(100.0f64)
    }
}
