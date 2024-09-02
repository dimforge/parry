use super::EPS;
use crate::math::{Point, Real, Vector};
use crate::query;
use crate::shape::{FeatureId, Segment, Triangle};
use crate::transformation::polygon_intersection::PolylinePointLocation;
use crate::utils::WBasis;

#[derive(Copy, Clone, Debug, Default)]
pub struct TriangleTriangleIntersectionPoint {
    pub p1: Point<Real>,
}

#[derive(Clone, Debug)]
pub enum TriangleTriangleIntersection {
    Segment {
        a: TriangleTriangleIntersectionPoint,
        b: TriangleTriangleIntersectionPoint,
    },
    Polygon(Vec<TriangleTriangleIntersectionPoint>),
}

impl Default for TriangleTriangleIntersection {
    fn default() -> Self {
        Self::Segment {
            a: Default::default(),
            b: Default::default(),
        }
    }
}

pub fn triangle_triangle_intersection(
    tri1: &Triangle,
    tri2: &Triangle,
) -> Option<TriangleTriangleIntersection> {
    let normal1 = tri1.robust_normal();
    let normal2 = tri2.robust_normal();

    if let Some(intersection_dir) = normal1.cross(&normal2).try_normalize(1.0e-6) {
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

        for hit1 in hits1.into_iter().flatten() {
            if hit1.0 < range1[0].0 {
                range1[0] = hit1;
            }
            if hit1.0 > range1[1].0 {
                range1[1] = hit1;
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

        for hit2 in hits2.into_iter().flatten() {
            if hit2.0 < range2[0].0 {
                range2[0] = hit2;
            }
            if hit2.0 > range2[1].0 {
                range2[1] = hit2;
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

        let a = if range2[0].0 > range1[0].0 + EPS {
            TriangleTriangleIntersectionPoint { p1: range2[0].1 }
        } else {
            TriangleTriangleIntersectionPoint { p1: range1[0].1 }
        };

        let b = if range2[1].0 < range1[1].0 - EPS {
            TriangleTriangleIntersectionPoint { p1: range2[1].1 }
        } else {
            TriangleTriangleIntersectionPoint { p1: range1[1].1 }
        };

        Some(TriangleTriangleIntersection::Segment { a, b })
    } else {
        let unit_normal2 = normal2.normalize();
        if (tri1.a - tri2.a).dot(&unit_normal2) < EPS {
            let basis = unit_normal2.orthonormal_basis();
            let proj =
                |vect: Vector<Real>| na::Point2::new(vect.dot(&basis[0]), vect.dot(&basis[1]));

            let mut intersections = vec![];

            let pts1 = tri1.vertices();
            let pts2 = tri2.vertices();
            let poly1 = [
                proj(tri1.a - tri2.a),
                proj(tri1.b - tri2.a),
                proj(tri1.c - tri2.a),
            ];
            let poly2 = [
                proj(Vector::zeros()), // = proj(tri2.a - tri2.a)
                proj(tri2.b - tri2.a),
                proj(tri2.c - tri2.a),
            ];

            let convert_loc = |loc, pts: &[Point<Real>; 3]| match loc {
                PolylinePointLocation::OnVertex(vid) => (FeatureId::Vertex(vid as u32), pts[vid]),
                PolylinePointLocation::OnEdge(vid1, vid2, bcoords) => (
                    match (vid1, vid2) {
                        (0, 1) | (1, 0) => FeatureId::Edge(0),
                        (1, 2) | (2, 1) => FeatureId::Edge(1),
                        (2, 0) | (0, 2) => FeatureId::Edge(2),
                        _ => unreachable!(),
                    },
                    pts[vid1] * bcoords[0] + pts[vid2].coords * bcoords[1],
                ),
            };

            crate::transformation::convex_polygons_intersection(&poly1, &poly2, |pt1, pt2| {
                let intersection = match (pt1, pt2) {
                    (Some(loc1), Some(loc2)) => {
                        let (_f1, p1) = convert_loc(loc1, pts1);
                        let (_f2, _p2) = convert_loc(loc2, pts2);
                        TriangleTriangleIntersectionPoint { p1 }
                    }
                    (Some(loc1), None) => {
                        let (_f1, p1) = convert_loc(loc1, pts1);
                        TriangleTriangleIntersectionPoint { p1 }
                    }
                    (None, Some(loc2)) => {
                        let (_f2, p2) = convert_loc(loc2, pts2);
                        TriangleTriangleIntersectionPoint { p1: p2 }
                    }
                    (None, None) => unreachable!(),
                };
                intersections.push(intersection);
            });

            Some(TriangleTriangleIntersection::Polygon(intersections))
        } else {
            None
        }
    }
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

    let time_of_impact =
        query::details::line_toi_with_halfspace(plane_center, plane_normal, &segment.a, &dir)?;
    let scaled_toi = time_of_impact * dir_norm;

    if scaled_toi < -EPS || scaled_toi > dir_norm + EPS {
        None
    } else if scaled_toi <= EPS {
        Some((segment.a, FeatureId::Vertex(vids.0)))
    } else if scaled_toi >= dir_norm - EPS {
        Some((segment.b, FeatureId::Vertex(vids.1)))
    } else {
        Some((segment.a + dir * time_of_impact, FeatureId::Edge(eid)))
    }
}
