use super::EPS;
use crate::math::{Point, Real, Vector};
use crate::query;
use crate::shape::{FeatureId, Segment, Triangle};

#[derive(Copy, Clone, Debug, Default)]
pub struct TriangleTriangleIntersectionPoint {
    pub pt: Point<Real>,
    pub f1: FeatureId,
    pub f2: FeatureId,
}

#[derive(Copy, Clone, Debug, Default)]
pub struct TriangleTriangleIntersection {
    pub a: TriangleTriangleIntersectionPoint,
    pub b: TriangleTriangleIntersectionPoint,
}

pub fn triangle_triangle_intersection(
    tri1: &Triangle,
    tri2: &Triangle,
) -> Option<TriangleTriangleIntersection> {
    let normal1 = tri1.scaled_normal();
    let normal2 = tri2.scaled_normal();

    let intersection_dir = normal1.cross(&normal2).normalize();

    let mut range1 = [
        (Real::MAX, Point::origin(), FeatureId::Unknown),
        (-Real::MAX, Point::origin(), FeatureId::Unknown),
    ];
    let mut range2 = [
        (Real::MAX, Point::origin(), FeatureId::Unknown),
        (-Real::MAX, Point::origin(), FeatureId::Unknown),
    ];

    let hits1 = [
        segment_plane_intersection(&tri2.a, &normal2, &Segment::new(tri1.a, tri1.b), 0, (0, 1))
            .map(|(p, feat)| (intersection_dir.dot(&p.coords), p, feat)),
        segment_plane_intersection(&tri2.a, &normal2, &Segment::new(tri1.b, tri1.c), 1, (1, 2))
            .map(|(p, feat)| (intersection_dir.dot(&p.coords), p, feat)),
        segment_plane_intersection(&tri2.a, &normal2, &Segment::new(tri1.c, tri1.a), 2, (2, 0))
            .map(|(p, feat)| (intersection_dir.dot(&p.coords), p, feat)),
    ];

    for k in 0..3 {
        if let Some(hit1) = &hits1[k] {
            if hit1.0 < range1[0].0 {
                range1[0] = *hit1;
            }
            if hit1.0 > range1[1].0 {
                range1[1] = *hit1;
            }
        }
    }

    if range1[0].0 >= range1[1].0 {
        // The first triangle doesn’t intersect the second plane.
        return None;
    }

    let hits2 = [
        segment_plane_intersection(&tri1.a, &normal1, &Segment::new(tri2.a, tri2.b), 0, (0, 1))
            .map(|(p, feat)| (intersection_dir.dot(&p.coords), p, feat)),
        segment_plane_intersection(&tri1.a, &normal1, &Segment::new(tri2.b, tri2.c), 1, (1, 2))
            .map(|(p, feat)| (intersection_dir.dot(&p.coords), p, feat)),
        segment_plane_intersection(&tri1.a, &normal1, &Segment::new(tri2.c, tri2.a), 2, (2, 0))
            .map(|(p, feat)| (intersection_dir.dot(&p.coords), p, feat)),
    ];

    for k in 0..3 {
        if let Some(hit2) = &hits2[k] {
            if hit2.0 < range2[0].0 {
                range2[0] = *hit2;
            }
            if hit2.0 > range2[1].0 {
                range2[1] = *hit2;
            }
        }
    }

    if range2[0].0 >= range2[1].0 {
        // The second triangle doesn’t intersect the first plane.
        return None;
    }

    if range1[1].0 <= range2[0].0 + EPS || range2[1].0 <= range1[0].0 + EPS {
        // The two triangles intersect each others’ plane, but these intersections are disjoint.
        return None;
    }

    let edge_between = |a, b| match (a, b) {
        (0, 1) | (1, 0) => FeatureId::Edge(0),
        (1, 2) | (2, 1) => FeatureId::Edge(1),
        (2, 0) | (0, 2) => FeatureId::Edge(2),
        _ => FeatureId::Edge(a),
    };

    let inter_f1 = match (range1[0].2, range1[1].2) {
        (FeatureId::Vertex(a), FeatureId::Vertex(b)) => edge_between(a, b),
        (FeatureId::Vertex(v), FeatureId::Edge(e)) | (FeatureId::Edge(e), FeatureId::Vertex(v)) => {
            if e == (v + 1) % 3 {
                FeatureId::Face(0)
            } else {
                FeatureId::Edge(e)
            }
        }
        _ => FeatureId::Face(0),
    };
    let inter_f2 = match (range2[0].2, range2[1].2) {
        (FeatureId::Vertex(a), FeatureId::Vertex(b)) => edge_between(a, b),
        (FeatureId::Vertex(v), FeatureId::Edge(e)) | (FeatureId::Edge(e), FeatureId::Vertex(v)) => {
            if e == (v + 1) % 3 {
                FeatureId::Face(0)
            } else {
                FeatureId::Edge(e)
            }
        }
        _ => FeatureId::Face(0),
    };

    let a = if range2[0].0 > range1[0].0 + EPS {
        TriangleTriangleIntersectionPoint {
            pt: range2[0].1,
            f1: inter_f1,
            f2: range2[0].2,
        }
    } else if range2[0].0 < range1[0].0 - EPS {
        TriangleTriangleIntersectionPoint {
            pt: range1[0].1,
            f1: range1[0].2,
            f2: inter_f2,
        }
    } else {
        TriangleTriangleIntersectionPoint {
            pt: range1[0].1,
            f1: range1[0].2,
            f2: range2[0].2,
        }
    };

    let b = if range2[1].0 < range1[1].0 - EPS {
        TriangleTriangleIntersectionPoint {
            pt: range2[1].1,
            f1: inter_f1,
            f2: range2[1].2,
        }
    } else if range2[1].0 > range1[1].0 + EPS {
        TriangleTriangleIntersectionPoint {
            pt: range1[1].1,
            f1: range1[1].2,
            f2: inter_f2,
        }
    } else {
        TriangleTriangleIntersectionPoint {
            pt: range1[1].1,
            f1: range1[1].2,
            f2: range2[1].2,
        }
    };

    Some(TriangleTriangleIntersection { a, b })
}

fn segment_plane_intersection(
    plane_center: &Point<Real>,
    plane_normal: &Vector<Real>,
    segment: &Segment,
    eid: u32,
    vids: (u32, u32),
) -> Option<(Point<Real>, FeatureId)> {
    let dir = segment.b - segment.a;
    let dir_norm = dir.norm();

    let toi =
        query::details::line_toi_with_halfspace(plane_center, &plane_normal, &segment.a, &dir)?;
    let scaled_toi = toi * dir_norm;

    if scaled_toi < -EPS || scaled_toi > dir_norm + EPS {
        None
    } else {
        if scaled_toi <= EPS {
            Some((segment.a, FeatureId::Vertex(vids.0)))
        } else if scaled_toi >= dir_norm - EPS {
            Some((segment.b, FeatureId::Vertex(vids.1)))
        } else {
            Some((segment.a + dir * toi, FeatureId::Edge(eid)))
        }
    }
}
