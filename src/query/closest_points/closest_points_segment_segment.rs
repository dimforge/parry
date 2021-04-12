use crate::math::{Isometry, Real};
use crate::query::ClosestPoints;
use crate::shape::{Segment, SegmentPointLocation};

use na::{self, Point};

/// Closest points between segments.
#[inline]
pub fn closest_points_segment_segment(
    pos12: &Isometry<Real>,
    seg1: &Segment,
    seg2: &Segment,
    margin: Real,
) -> ClosestPoints {
    let (loc1, loc2) = closest_points_segment_segment_with_locations(pos12, seg1, seg2);
    let p1 = seg1.point_at(&loc1);
    let p2 = seg2.point_at(&loc2);

    if na::distance_squared(&p1, &(pos12 * p2)) <= margin * margin {
        ClosestPoints::WithinMargin(p1, p2)
    } else {
        ClosestPoints::Disjoint
    }
}

// FIXME: use this specialized procedure for distance/interference/contact determination as well.
/// Closest points between two segments.
#[inline]
pub fn closest_points_segment_segment_with_locations(
    pos12: &Isometry<Real>,
    seg1: &Segment,
    seg2: &Segment,
) -> (SegmentPointLocation, SegmentPointLocation) {
    let seg2_1 = seg2.transformed(pos12);
    closest_points_segment_segment_with_locations_nD((&seg1.a, &seg1.b), (&seg2_1.a, &seg2_1.b))
}

/// Segment-segment closest points computation in an arbitrary dimension.
#[allow(non_snake_case)]
#[inline]
pub fn closest_points_segment_segment_with_locations_nD<const D: usize>(
    seg1: (&Point<Real, D>, &Point<Real, D>),
    seg2: (&Point<Real, D>, &Point<Real, D>),
) -> (SegmentPointLocation, SegmentPointLocation) {
    // Inspired by RealField-time collision detection by Christer Ericson.
    let d1 = seg1.1 - seg1.0;
    let d2 = seg2.1 - seg2.0;
    let r = seg1.0 - seg2.0;

    let a = d1.norm_squared();
    let e = d2.norm_squared();
    let f = d2.dot(&r);

    let _0: Real = 0.0;
    let _1: Real = 1.0;

    let mut s;
    let mut t;

    let _eps = crate::math::DEFAULT_EPSILON;
    if a <= _eps && e <= _eps {
        s = _0;
        t = _0;
    } else if a <= _eps {
        s = _0;
        t = na::clamp(f / e, _0, _1);
    } else {
        let c = d1.dot(&r);
        if e <= _eps {
            t = _0;
            s = na::clamp(-c / a, _0, _1);
        } else {
            let b = d1.dot(&d2);
            let ae = a * e;
            let bb = b * b;
            let denom = ae - bb;

            // Use absolute and ulps error to test collinearity.
            if denom > _eps && !ulps_eq!(ae, bb) {
                s = na::clamp((b * f - c * e) / denom, _0, _1);
            } else {
                s = _0;
            }

            t = (b * s + f) / e;

            if t < _0 {
                t = _0;
                s = na::clamp(-c / a, _0, _1);
            } else if t > _1 {
                t = _1;
                s = na::clamp((b - c) / a, _0, _1);
            }
        }
    }

    let loc1 = if s == _0 {
        SegmentPointLocation::OnVertex(0)
    } else if s == _1 {
        SegmentPointLocation::OnVertex(1)
    } else {
        SegmentPointLocation::OnEdge([_1 - s, s])
    };

    let loc2 = if t == _0 {
        SegmentPointLocation::OnVertex(0)
    } else if t == _1 {
        SegmentPointLocation::OnVertex(1)
    } else {
        SegmentPointLocation::OnEdge([_1 - t, t])
    };

    (loc1, loc2)
}
