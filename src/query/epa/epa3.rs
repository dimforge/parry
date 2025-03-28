//! Three-dimensional penetration depth queries using the Expanding Polytope Algorithm.

use crate::math::{Isometry, Point, Real, Vector};
use crate::query::gjk::{self, CSOPoint, ConstantOrigin, VoronoiSimplex};
use crate::query::PointQueryWithLocation;
use crate::shape::{SupportMap, Triangle, TrianglePointLocation};
use crate::utils;
use alloc::collections::BinaryHeap;
use alloc::vec::Vec;
use core::cmp::Ordering;
use na::{self, Unit};
use num::Bounded;

#[derive(Copy, Clone, PartialEq)]
struct FaceId {
    id: usize,
    neg_dist: Real,
}

impl FaceId {
    fn new(id: usize, neg_dist: Real) -> Option<Self> {
        if neg_dist > gjk::eps_tol() {
            None
        } else {
            Some(FaceId { id, neg_dist })
        }
    }
}

impl Eq for FaceId {}

impl PartialOrd for FaceId {
    #[inline]
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for FaceId {
    #[inline]
    fn cmp(&self, other: &Self) -> Ordering {
        if self.neg_dist < other.neg_dist {
            Ordering::Less
        } else if self.neg_dist > other.neg_dist {
            Ordering::Greater
        } else {
            Ordering::Equal
        }
    }
}

#[derive(Clone, Debug)]
struct Face {
    pts: [usize; 3],
    adj: [usize; 3],
    normal: Unit<Vector<Real>>,
    bcoords: [Real; 3],
    deleted: bool,
}

impl Face {
    pub fn new_with_proj(
        vertices: &[CSOPoint],
        bcoords: [Real; 3],
        pts: [usize; 3],
        adj: [usize; 3],
    ) -> Self {
        let normal;

        if let Some(n) = utils::ccw_face_normal([
            &vertices[pts[0]].point,
            &vertices[pts[1]].point,
            &vertices[pts[2]].point,
        ]) {
            normal = n;
        } else {
            // This is a bit of a hack for degenerate faces.
            // TODO: It will work OK with our current code, though
            // we should do this in another way to avoid any risk
            // of misusing the face normal in the future.
            normal = Unit::new_unchecked(na::zero());
        }

        Face {
            pts,
            bcoords,
            adj,
            normal,
            deleted: false,
        }
    }

    pub fn new(vertices: &[CSOPoint], pts: [usize; 3], adj: [usize; 3]) -> (Self, bool) {
        let tri = Triangle::new(
            vertices[pts[0]].point,
            vertices[pts[1]].point,
            vertices[pts[2]].point,
        );
        let (proj, loc) = tri.project_local_point_and_get_location(&Point::<Real>::origin(), true);

        match loc {
            TrianglePointLocation::OnVertex(_) | TrianglePointLocation::OnEdge(_, _) => {
                let eps_tol = crate::math::DEFAULT_EPSILON * 100.0; // Same as in closest_points
                (
                    // barycentric_coordinates is guaranteed to work in OnVertex and OnEdge locations
                    Self::new_with_proj(vertices, loc.barycentric_coordinates().unwrap(), pts, adj),
                    proj.is_inside_eps(&Point::<Real>::origin(), eps_tol),
                )
            }
            TrianglePointLocation::OnFace(_, bcoords) => {
                (Self::new_with_proj(vertices, bcoords, pts, adj), true)
            }
            _ => (Self::new_with_proj(vertices, [0.0; 3], pts, adj), false),
        }
    }

    pub fn closest_points(&self, vertices: &[CSOPoint]) -> (Point<Real>, Point<Real>) {
        (
            vertices[self.pts[0]].orig1 * self.bcoords[0]
                + vertices[self.pts[1]].orig1.coords * self.bcoords[1]
                + vertices[self.pts[2]].orig1.coords * self.bcoords[2],
            vertices[self.pts[0]].orig2 * self.bcoords[0]
                + vertices[self.pts[1]].orig2.coords * self.bcoords[1]
                + vertices[self.pts[2]].orig2.coords * self.bcoords[2],
        )
    }

    pub fn contains_point(&self, id: usize) -> bool {
        self.pts[0] == id || self.pts[1] == id || self.pts[2] == id
    }

    pub fn next_ccw_pt_id(&self, id: usize) -> usize {
        if self.pts[0] == id {
            1
        } else if self.pts[1] == id {
            2
        } else {
            if self.pts[2] != id {
                log::debug!(
                    "Hit unexpected state in EPA: found index {}, expected: {}.",
                    self.pts[2],
                    id
                );
            }

            0
        }
    }

    pub fn can_be_seen_by(&self, vertices: &[CSOPoint], point: usize, opp_pt_id: usize) -> bool {
        let p0 = &vertices[self.pts[opp_pt_id]].point;
        let p1 = &vertices[self.pts[(opp_pt_id + 1) % 3]].point;
        let p2 = &vertices[self.pts[(opp_pt_id + 2) % 3]].point;
        let pt = &vertices[point].point;

        // NOTE: it is important that we return true for the case where
        // the dot product is zero. This is because degenerate faces will
        // have a zero normal, causing the dot product to be zero.
        // So return true for these case will let us skip the triangle
        // during silhouette computation.
        (*pt - *p0).dot(&self.normal) >= -gjk::eps_tol()
            || Triangle::new(*p1, *p2, *pt).is_affinely_dependent()
    }
}

struct SilhouetteEdge {
    face_id: usize,
    opp_pt_id: usize,
}

impl SilhouetteEdge {
    pub fn new(face_id: usize, opp_pt_id: usize) -> Self {
        SilhouetteEdge { face_id, opp_pt_id }
    }
}

/// The Expanding Polytope Algorithm in 3D.
#[derive(Default)]
pub struct EPA {
    vertices: Vec<CSOPoint>,
    faces: Vec<Face>,
    silhouette: Vec<SilhouetteEdge>,
    heap: BinaryHeap<FaceId>,
}

impl EPA {
    /// Creates a new instance of the 3D Expanding Polytope Algorithm.
    pub fn new() -> Self {
        Self::default()
    }

    fn reset(&mut self) {
        self.vertices.clear();
        self.faces.clear();
        self.heap.clear();
        self.silhouette.clear();
    }

    /// Projects the origin on boundary of the given shape.
    ///
    /// The origin is assumed to be inside of the shape. If it is outside
    /// use the GJK algorithm instead.
    /// Return `None` if the origin is not inside of the shape or if
    /// the EPA algorithm failed to compute the projection.
    ///
    /// Return the projected point in the local-space of `g`.
    pub fn project_origin<G: ?Sized + SupportMap>(
        &mut self,
        m: &Isometry<Real>,
        g: &G,
        simplex: &VoronoiSimplex,
    ) -> Option<Point<Real>> {
        self.closest_points(&m.inverse(), g, &ConstantOrigin, simplex)
            .map(|(p, _, _)| p)
    }

    /// Projects the origin on a shape using the EPA algorithm.
    ///
    /// The origin is assumed to be located inside of the shape.
    /// Returns `None` if the EPA fails to converge or if `g1` and `g2` are not penetrating.
    pub fn closest_points<G1, G2>(
        &mut self,
        pos12: &Isometry<Real>,
        g1: &G1,
        g2: &G2,
        simplex: &VoronoiSimplex,
    ) -> Option<(Point<Real>, Point<Real>, Unit<Vector<Real>>)>
    where
        G1: ?Sized + SupportMap,
        G2: ?Sized + SupportMap,
    {
        let _eps = crate::math::DEFAULT_EPSILON;
        let _eps_tol = _eps * 100.0;

        self.reset();

        /*
         * Initialization.
         */
        for i in 0..simplex.dimension() + 1 {
            self.vertices.push(*simplex.point(i));
        }

        if simplex.dimension() == 0 {
            let mut n: Vector<Real> = na::zero();
            n[1] = 1.0;
            return Some((Point::origin(), Point::origin(), Unit::new_unchecked(n)));
        } else if simplex.dimension() == 3 {
            let dp1 = self.vertices[1] - self.vertices[0];
            let dp2 = self.vertices[2] - self.vertices[0];
            let dp3 = self.vertices[3] - self.vertices[0];

            if dp1.cross(&dp2).dot(&dp3) > 0.0 {
                self.vertices.swap(1, 2)
            }

            let pts1 = [0, 1, 2];
            let pts2 = [1, 3, 2];
            let pts3 = [0, 2, 3];
            let pts4 = [0, 3, 1];

            let adj1 = [3, 1, 2];
            let adj2 = [3, 2, 0];
            let adj3 = [0, 1, 3];
            let adj4 = [2, 1, 0];

            let (face1, proj_inside1) = Face::new(&self.vertices, pts1, adj1);
            let (face2, proj_inside2) = Face::new(&self.vertices, pts2, adj2);
            let (face3, proj_inside3) = Face::new(&self.vertices, pts3, adj3);
            let (face4, proj_inside4) = Face::new(&self.vertices, pts4, adj4);

            self.faces.push(face1);
            self.faces.push(face2);
            self.faces.push(face3);
            self.faces.push(face4);

            if proj_inside1 {
                let dist1 = self.faces[0].normal.dot(&self.vertices[0].point.coords);
                self.heap.push(FaceId::new(0, -dist1)?);
            }

            if proj_inside2 {
                let dist2 = self.faces[1].normal.dot(&self.vertices[1].point.coords);
                self.heap.push(FaceId::new(1, -dist2)?);
            }

            if proj_inside3 {
                let dist3 = self.faces[2].normal.dot(&self.vertices[2].point.coords);
                self.heap.push(FaceId::new(2, -dist3)?);
            }

            if proj_inside4 {
                let dist4 = self.faces[3].normal.dot(&self.vertices[3].point.coords);
                self.heap.push(FaceId::new(3, -dist4)?);
            }

            if !(proj_inside1 || proj_inside2 || proj_inside3 || proj_inside4) {
                // Related issues:
                // https://github.com/dimforge/parry/issues/253
                // https://github.com/dimforge/parry/issues/246
                log::debug!("Hit unexpected state in EPA: failed to project the origin on the initial simplex.");
                return None;
            }
        } else {
            if simplex.dimension() == 1 {
                let dpt = self.vertices[1] - self.vertices[0];

                Vector::orthonormal_subspace_basis(&[dpt], |dir| {
                    // XXX: dir should already be unit on nalgebra!
                    let dir = Unit::new_unchecked(*dir);
                    self.vertices
                        .push(CSOPoint::from_shapes(pos12, g1, g2, &dir));
                    false
                });
            }

            let pts1 = [0, 1, 2];
            let pts2 = [0, 2, 1];

            let adj1 = [1, 1, 1];
            let adj2 = [0, 0, 0];

            let (face1, _) = Face::new(&self.vertices, pts1, adj1);
            let (face2, _) = Face::new(&self.vertices, pts2, adj2);
            self.faces.push(face1);
            self.faces.push(face2);

            self.heap.push(FaceId::new(0, 0.0)?);
            self.heap.push(FaceId::new(1, 0.0)?);
        }

        let mut niter = 0;
        let mut max_dist = Real::max_value();
        let mut best_face_id = *self.heap.peek()?;
        let mut old_dist = 0.0;

        /*
         * Run the expansion.
         */
        while let Some(face_id) = self.heap.pop() {
            // Create new faces.
            let face = self.faces[face_id.id].clone();

            if face.deleted {
                continue;
            }

            let cso_point = CSOPoint::from_shapes(pos12, g1, g2, &face.normal);
            let support_point_id = self.vertices.len();
            self.vertices.push(cso_point);

            let candidate_max_dist = cso_point.point.coords.dot(&face.normal);

            if candidate_max_dist < max_dist {
                best_face_id = face_id;
                max_dist = candidate_max_dist;
            }

            let curr_dist = -face_id.neg_dist;

            if max_dist - curr_dist < _eps_tol ||
                // Accept the intersection as the algorithm is stuck and no new points will be found
                // This happens because of numerical stability issue
                ((curr_dist - old_dist).abs() < _eps && candidate_max_dist < max_dist)
            {
                let best_face = &self.faces[best_face_id.id];
                let points = best_face.closest_points(&self.vertices);
                return Some((points.0, points.1, best_face.normal));
            }

            old_dist = curr_dist;

            self.faces[face_id.id].deleted = true;

            let adj_opp_pt_id1 = self.faces[face.adj[0]].next_ccw_pt_id(face.pts[0]);
            let adj_opp_pt_id2 = self.faces[face.adj[1]].next_ccw_pt_id(face.pts[1]);
            let adj_opp_pt_id3 = self.faces[face.adj[2]].next_ccw_pt_id(face.pts[2]);

            self.compute_silhouette(support_point_id, face.adj[0], adj_opp_pt_id1);
            self.compute_silhouette(support_point_id, face.adj[1], adj_opp_pt_id2);
            self.compute_silhouette(support_point_id, face.adj[2], adj_opp_pt_id3);

            let first_new_face_id = self.faces.len();

            if self.silhouette.is_empty() {
                // TODO: Something went very wrong because we failed to extract a silhouette…
                return None;
            }

            for edge in &self.silhouette {
                if !self.faces[edge.face_id].deleted {
                    let new_face_id = self.faces.len();

                    let face_adj = &mut self.faces[edge.face_id];
                    let pt_id1 = face_adj.pts[(edge.opp_pt_id + 2) % 3];
                    let pt_id2 = face_adj.pts[(edge.opp_pt_id + 1) % 3];

                    let pts = [pt_id1, pt_id2, support_point_id];
                    let adj = [edge.face_id, new_face_id + 1, new_face_id - 1];
                    let new_face = Face::new(&self.vertices, pts, adj);

                    face_adj.adj[(edge.opp_pt_id + 1) % 3] = new_face_id;

                    self.faces.push(new_face.0);

                    if new_face.1 {
                        let pt = self.vertices[self.faces[new_face_id].pts[0]].point.coords;
                        let dist = self.faces[new_face_id].normal.dot(&pt);
                        if dist < curr_dist {
                            // TODO: if we reach this point, there were issues due to
                            // numerical errors.
                            let points = face.closest_points(&self.vertices);
                            return Some((points.0, points.1, face.normal));
                        }

                        self.heap.push(FaceId::new(new_face_id, -dist)?);
                    }
                }
            }

            if first_new_face_id == self.faces.len() {
                // Something went very wrong because all the edges
                // from the silhouette belonged to deleted faces.
                return None;
            }

            self.faces[first_new_face_id].adj[2] = self.faces.len() - 1;
            self.faces.last_mut().unwrap().adj[1] = first_new_face_id;

            self.silhouette.clear();
            // self.check_topology(); // NOTE: for debugging only.

            niter += 1;
            if niter > 100 {
                // if we reached this point, our algorithm didn't converge to what precision we wanted.
                // still return an intersection point, as it's probably close enough.
                break;
            }
        }

        let best_face = &self.faces[best_face_id.id];
        let points = best_face.closest_points(&self.vertices);
        Some((points.0, points.1, best_face.normal))
    }

    fn compute_silhouette(&mut self, point: usize, id: usize, opp_pt_id: usize) {
        if !self.faces[id].deleted {
            if !self.faces[id].can_be_seen_by(&self.vertices, point, opp_pt_id) {
                self.silhouette.push(SilhouetteEdge::new(id, opp_pt_id));
            } else {
                self.faces[id].deleted = true;

                let adj_pt_id1 = (opp_pt_id + 2) % 3;
                let adj_pt_id2 = opp_pt_id;

                let adj1 = self.faces[id].adj[adj_pt_id1];
                let adj2 = self.faces[id].adj[adj_pt_id2];

                let adj_opp_pt_id1 =
                    self.faces[adj1].next_ccw_pt_id(self.faces[id].pts[adj_pt_id1]);
                let adj_opp_pt_id2 =
                    self.faces[adj2].next_ccw_pt_id(self.faces[id].pts[adj_pt_id2]);

                self.compute_silhouette(point, adj1, adj_opp_pt_id1);
                self.compute_silhouette(point, adj2, adj_opp_pt_id2);
            }
        }
    }

    #[allow(dead_code)]
    #[cfg(feature = "std")]
    fn print_silhouette(&self) {
        use std::{print, println};

        print!("Silhouette points: ");
        for i in 0..self.silhouette.len() {
            let edge = &self.silhouette[i];
            let face = &self.faces[edge.face_id];

            if !face.deleted {
                print!(
                    "({}, {}) ",
                    face.pts[(edge.opp_pt_id + 2) % 3],
                    face.pts[(edge.opp_pt_id + 1) % 3]
                );
            }
        }
        println!();
    }

    #[allow(dead_code)]
    fn check_topology(&self) {
        for i in 0..self.faces.len() {
            let face = &self.faces[i];
            if face.deleted {
                continue;
            }

            // println!("checking {}-th face.", i);
            let adj1 = &self.faces[face.adj[0]];
            let adj2 = &self.faces[face.adj[1]];
            let adj3 = &self.faces[face.adj[2]];

            assert!(!adj1.deleted);
            assert!(!adj2.deleted);
            assert!(!adj3.deleted);

            assert!(face.pts[0] != face.pts[1]);
            assert!(face.pts[0] != face.pts[2]);
            assert!(face.pts[1] != face.pts[2]);

            assert!(adj1.contains_point(face.pts[0]));
            assert!(adj1.contains_point(face.pts[1]));

            assert!(adj2.contains_point(face.pts[1]));
            assert!(adj2.contains_point(face.pts[2]));

            assert!(adj3.contains_point(face.pts[2]));
            assert!(adj3.contains_point(face.pts[0]));

            let opp_pt_id1 = adj1.next_ccw_pt_id(face.pts[0]);
            let opp_pt_id2 = adj2.next_ccw_pt_id(face.pts[1]);
            let opp_pt_id3 = adj3.next_ccw_pt_id(face.pts[2]);

            assert!(!face.contains_point(adj1.pts[opp_pt_id1]));
            assert!(!face.contains_point(adj2.pts[opp_pt_id2]));
            assert!(!face.contains_point(adj3.pts[opp_pt_id3]));

            assert!(adj1.adj[(opp_pt_id1 + 1) % 3] == i);
            assert!(adj2.adj[(opp_pt_id2 + 1) % 3] == i);
            assert!(adj3.adj[(opp_pt_id3 + 1) % 3] == i);
        }
    }
}
