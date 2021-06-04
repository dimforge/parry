//! Definition of the triangle shape.

use crate::math::{Isometry, Point, Real, Vector};
use crate::shape::{FeatureId, SupportMap};
use crate::shape::{PolygonalFeature, Segment};
use crate::utils;

use na::{self, ComplexField, Unit};
use num::Zero;
#[cfg(feature = "dim3")]
use std::f64;
use std::mem;

/// A triangle shape.
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[repr(C)]
#[derive(PartialEq, Debug, Copy, Clone)]
pub struct Triangle {
    /// The triangle first point.
    pub a: Point<Real>,
    /// The triangle second point.
    pub b: Point<Real>,
    /// The triangle third point.
    pub c: Point<Real>,
}

/// Description of the location of a point on a triangle.
#[derive(Copy, Clone, Debug)]
pub enum TrianglePointLocation {
    /// The point lies on a vertex.
    OnVertex(u32),
    /// The point lies on an edge.
    ///
    /// The 0-st edge is the segment AB.
    /// The 1-st edge is the segment BC.
    /// The 2-nd edge is the segment AC.
    // XXX: it appears the conversion of edge indexing here does not match the
    // convension of edge indexing for the `fn edge` method (from the ConvexPolyhedron impl).
    OnEdge(u32, [Real; 2]),
    /// The point lies on the triangle interior.
    ///
    /// The integer indicates on which side of the face the point is. 0 indicates the point
    /// is on the half-space toward the CW normal of the triangle. 1 indicates the point is on the other
    /// half-space. This is always set to 0 in 2D.
    OnFace(u32, [Real; 3]),
    /// The point lies on the triangle interior (for "solid" point queries).
    OnSolid,
}

impl TrianglePointLocation {
    /// The barycentric coordinates corresponding to this point location.
    ///
    /// Returns `None` if the location is `TrianglePointLocation::OnSolid`.
    pub fn barycentric_coordinates(&self) -> Option<[Real; 3]> {
        let mut bcoords = [0.0; 3];

        match self {
            TrianglePointLocation::OnVertex(i) => bcoords[*i as usize] = 1.0,
            TrianglePointLocation::OnEdge(i, uv) => {
                let idx = match i {
                    0 => (0, 1),
                    1 => (1, 2),
                    2 => (0, 2),
                    _ => unreachable!(),
                };

                bcoords[idx.0] = uv[0];
                bcoords[idx.1] = uv[1];
            }
            TrianglePointLocation::OnFace(_, uvw) => {
                bcoords[0] = uvw[0];
                bcoords[1] = uvw[1];
                bcoords[2] = uvw[2];
            }
            TrianglePointLocation::OnSolid => {
                return None;
            }
        }

        Some(bcoords)
    }

    /// Returns `true` if the point is located on the relative interior of the triangle.
    pub fn is_on_face(&self) -> bool {
        if let TrianglePointLocation::OnFace(..) = *self {
            true
        } else {
            false
        }
    }
}

impl From<[Point<Real>; 3]> for Triangle {
    fn from(arr: [Point<Real>; 3]) -> Self {
        *Self::from_array(&arr)
    }
}

impl Triangle {
    /// Creates a triangle from three points.
    #[inline]
    pub fn new(a: Point<Real>, b: Point<Real>, c: Point<Real>) -> Triangle {
        Triangle { a, b, c }
    }

    /// Creates the reference to a triangle from the reference to an array of three points.
    pub fn from_array(arr: &[Point<Real>; 3]) -> &Triangle {
        unsafe { mem::transmute(arr) }
    }

    /// Reference to an array containing the three vertices of this triangle.
    #[inline]
    pub fn vertices(&self) -> &[Point<Real>; 3] {
        unsafe { mem::transmute(self) }
    }

    /// The normal of this triangle assuming it is oriented ccw.
    ///
    /// The normal points such that it is collinear to `AB × AC` (where `×` denotes the cross
    /// product).
    #[inline]
    pub fn normal(&self) -> Option<Unit<Vector<Real>>> {
        Unit::try_new(self.scaled_normal(), crate::math::DEFAULT_EPSILON)
    }

    /// The three edges of this triangle: [AB, BC, CA].
    #[inline]
    pub fn edges(&self) -> [Segment; 3] {
        [
            Segment::new(self.a, self.b),
            Segment::new(self.b, self.c),
            Segment::new(self.c, self.a),
        ]
    }

    /// Returns a new triangle with vertices transformed by `m`.
    #[inline]
    pub fn transformed(&self, m: &Isometry<Real>) -> Self {
        Triangle::new(m * self.a, m * self.b, m * self.c)
    }
    /// The three edges scaled directions of this triangle: [B - A, C - B, A - C].
    #[inline]
    pub fn edges_scaled_directions(&self) -> [Vector<Real>; 3] {
        [self.b - self.a, self.c - self.b, self.a - self.c]
    }

    /// Return the edge segment of this cuboid with a normal cone containing
    /// a direction that that maximizes the dot product with `local_dir`.
    pub fn local_support_edge_segment(&self, dir: Vector<Real>) -> Segment {
        let dots = na::Vector3::new(
            dir.dot(&self.a.coords),
            dir.dot(&self.b.coords),
            dir.dot(&self.c.coords),
        );

        match dots.imin() {
            0 => Segment::new(self.b, self.c),
            1 => Segment::new(self.c, self.a),
            _ => Segment::new(self.a, self.b),
        }
    }

    /// Return the face of this triangle with a normal that maximizes
    /// the dot product with `dir`.
    #[cfg(feature = "dim3")]
    pub fn support_face(&self, _dir: Vector<Real>) -> PolygonalFeature {
        PolygonalFeature::from(*self)
    }

    /// Return the face of this triangle with a normal that maximizes
    /// the dot product with `dir`.
    #[cfg(feature = "dim2")]
    pub fn support_face(&self, dir: Vector<Real>) -> PolygonalFeature {
        let mut best = 0;
        let mut best_dot = -Real::MAX;

        for (i, tangent) in self.edges_scaled_directions().iter().enumerate() {
            let normal = Vector::new(tangent.y, -tangent.x);
            if let Some(normal) = Unit::try_new(normal, 0.0) {
                let dot = normal.dot(&dir);
                if normal.dot(&dir) > best_dot {
                    best = i;
                    best_dot = dot;
                }
            }
        }

        let pts = self.vertices();
        let i1 = best;
        let i2 = (best + 1) % 3;

        PolygonalFeature {
            vertices: [pts[i1], pts[i2]],
            vids: [i1 as u32, i2 as u32],
            fid: i1 as u32,
            num_vertices: 2,
        }
    }

    /// A vector normal of this triangle.
    ///
    /// The vector points such that it is collinear to `AB × AC` (where `×` denotes the cross
    /// product).
    #[inline]
    pub fn scaled_normal(&self) -> Vector<Real> {
        let ab = self.b - self.a;
        let ac = self.c - self.a;
        ab.cross(&ac)
    }

    /// Computes the extents of this triangle on the given direction.
    ///
    /// This computes the min and max values of the dot products between each
    /// vertex of this triangle and `dir`.
    #[inline]
    pub fn extents_on_dir(&self, dir: &Unit<Vector<Real>>) -> (Real, Real) {
        let a = self.a.coords.dot(dir);
        let b = self.b.coords.dot(dir);
        let c = self.c.coords.dot(dir);

        if a > b {
            if b > c {
                (c, a)
            } else if a > c {
                (b, a)
            } else {
                (b, c)
            }
        } else {
            // b >= a
            if a > c {
                (c, b)
            } else if b > c {
                (a, b)
            } else {
                (a, c)
            }
        }
    }
    //
    // #[cfg(feature = "dim3")]
    // fn support_feature_id_toward(&self, local_dir: &Unit<Vector<Real>>, eps: Real) -> FeatureId {
    //     if let Some(normal) = self.normal() {
    //         let (seps, ceps) = ComplexField::sin_cos(eps);
    //
    //         let normal_dot = local_dir.dot(&*normal);
    //         if normal_dot >= ceps {
    //             FeatureId::Face(0)
    //         } else if normal_dot <= -ceps {
    //             FeatureId::Face(1)
    //         } else {
    //             let edges = self.edges();
    //             let mut dots = [0.0; 3];
    //
    //             let dir1 = edges[0].direction();
    //             if let Some(dir1) = dir1 {
    //                 dots[0] = dir1.dot(local_dir);
    //
    //                 if dots[0].abs() < seps {
    //                     return FeatureId::Edge(0);
    //                 }
    //             }
    //
    //             let dir2 = edges[1].direction();
    //             if let Some(dir2) = dir2 {
    //                 dots[1] = dir2.dot(local_dir);
    //
    //                 if dots[1].abs() < seps {
    //                     return FeatureId::Edge(1);
    //                 }
    //             }
    //
    //             let dir3 = edges[2].direction();
    //             if let Some(dir3) = dir3 {
    //                 dots[2] = dir3.dot(local_dir);
    //
    //                 if dots[2].abs() < seps {
    //                     return FeatureId::Edge(2);
    //                 }
    //             }
    //
    //             if dots[0] > 0.0 && dots[1] < 0.0 {
    //                 FeatureId::Vertex(1)
    //             } else if dots[1] > 0.0 && dots[2] < 0.0 {
    //                 FeatureId::Vertex(2)
    //             } else {
    //                 FeatureId::Vertex(0)
    //             }
    //         }
    //     } else {
    //         FeatureId::Vertex(0)
    //     }
    // }

    /// The area of this triangle.
    #[inline]
    pub fn area(&self) -> Real {
        // Kahan's formula.
        let mut a = na::distance(&self.a, &self.b);
        let mut b = na::distance(&self.b, &self.c);
        let mut c = na::distance(&self.c, &self.a);

        let (c, b, a) = utils::sort3(&mut a, &mut b, &mut c);
        let a = *a;
        let b = *b;
        let c = *c;

        let sqr = (a + (b + c)) * (c - (a - b)) * (c + (a - b)) * (a + (b - c));

        // We take the max(0.0) because it can be slightly negative
        // because of numerical errors due to almost-degenerate triangles.
        ComplexField::sqrt(sqr.max(0.0)) * 0.25
    }

    /// Computes the unit angular inertia of this triangle.
    #[cfg(feature = "dim2")]
    pub fn unit_angular_inertia(&self) -> Real {
        let factor = 1.0 / 6.0;

        // Algorithm adapted from Box2D
        let e1 = self.b - self.a;
        let e2 = self.c - self.a;

        let intx2 = e1.x * e1.x + e2.x * e1.x + e2.x * e2.x;
        let inty2 = e1.y * e1.y + e2.y * e1.y + e2.y * e2.y;
        factor * (intx2 + inty2)
    }

    /// The geometric center of this triangle.
    #[inline]
    pub fn center(&self) -> Point<Real> {
        utils::center(&[self.a, self.b, self.c])
    }

    /// The perimeter of this triangle.
    #[inline]
    pub fn perimeter(&self) -> Real {
        na::distance(&self.a, &self.b)
            + na::distance(&self.b, &self.c)
            + na::distance(&self.c, &self.a)
    }

    /// The circumcircle of this triangle.
    pub fn circumcircle(&self) -> (Point<Real>, Real) {
        let a = self.a - self.c;
        let b = self.b - self.c;

        let na = a.norm_squared();
        let nb = b.norm_squared();

        let dab = a.dot(&b);

        let _2: Real = na::convert::<f64, Real>(2.0);
        let denom = _2 * (na * nb - dab * dab);

        if denom.is_zero() {
            // The triangle is degenerate (the three points are colinear).
            // So we find the longest segment and take its center.
            let c = self.a - self.b;
            let nc = c.norm_squared();

            if nc >= na && nc >= nb {
                // Longest segment: [&self.a, &self.b]
                (
                    na::center(&self.a, &self.b),
                    ComplexField::sqrt(nc) / na::convert::<f64, Real>(2.0f64),
                )
            } else if na >= nb && na >= nc {
                // Longest segment: [&self.a, pc]
                (
                    na::center(&self.a, &self.c),
                    ComplexField::sqrt(na) / na::convert::<f64, Real>(2.0f64),
                )
            } else {
                // Longest segment: [&self.b, &self.c]
                (
                    na::center(&self.b, &self.c),
                    ComplexField::sqrt(nb) / na::convert::<f64, Real>(2.0f64),
                )
            }
        } else {
            let k = b * na - a * nb;

            let center = self.c + (a * k.dot(&b) - b * k.dot(&a)) / denom;
            let radius = na::distance(&self.a, &center);

            (center, radius)
        }
    }

    /// Tests if this triangle is affinely dependent, i.e., its points are almost aligned.
    #[cfg(feature = "dim3")]
    pub fn is_affinely_dependent(&self) -> bool {
        const EPS: Real = crate::math::DEFAULT_EPSILON * 100.0;

        let p1p2 = self.b - self.a;
        let p1p3 = self.c - self.a;
        relative_eq!(p1p2.cross(&p1p3).norm_squared(), 0.0, epsilon = EPS * EPS)

        // relative_eq!(
        //     self.area(),
        //     0.0,
        //     epsilon = EPS * self.perimeter()
        // )
    }

    /// Tests if a point is inside of this triangle.
    pub fn contains_point(&self, p: &Point<Real>) -> bool {
        let p1p2 = self.b - self.a;
        let p2p3 = self.c - self.b;
        let p3p1 = self.a - self.c;

        let p1p = *p - self.a;
        let p2p = *p - self.b;
        let p3p = *p - self.c;

        let d11 = p1p.dot(&p1p2);
        let d12 = p2p.dot(&p2p3);
        let d13 = p3p.dot(&p3p1);

        d11 >= 0.0
            && d11 <= p1p2.norm_squared()
            && d12 >= 0.0
            && d12 <= p2p3.norm_squared()
            && d13 >= 0.0
            && d13 <= p3p1.norm_squared()
    }

    /// The normal of the given feature of this shape.
    pub fn feature_normal(&self, _: FeatureId) -> Option<Unit<Vector<Real>>> {
        self.normal()
    }
}

impl SupportMap for Triangle {
    #[inline]
    fn local_support_point(&self, dir: &Vector<Real>) -> Point<Real> {
        let d1 = self.a.coords.dot(dir);
        let d2 = self.b.coords.dot(dir);
        let d3 = self.c.coords.dot(dir);

        if d1 > d2 {
            if d1 > d3 {
                self.a
            } else {
                self.c
            }
        } else {
            if d2 > d3 {
                self.b
            } else {
                self.c
            }
        }
    }
}

/*
#[cfg(feature = "dim3")]
impl ConvexPolyhedron for Triangle {
    fn vertex(&self, id: FeatureId) -> Point<Real> {
        match id.unwrap_vertex() {
            0 => self.a,
            1 => self.b,
            2 => self.c,
            _ => panic!("Triangle vertex index out of bounds."),
        }
    }
    fn edge(&self, id: FeatureId) -> (Point<Real>, Point<Real>, FeatureId, FeatureId) {
        match id.unwrap_edge() {
            0 => (self.a, self.b, FeatureId::Vertex(0), FeatureId::Vertex(1)),
            1 => (self.b, self.c, FeatureId::Vertex(1), FeatureId::Vertex(2)),
            2 => (self.c, self.a, FeatureId::Vertex(2), FeatureId::Vertex(0)),
            _ => panic!("Triangle edge index out of bounds."),
        }
    }

    fn face(&self, id: FeatureId, face: &mut ConvexPolygonalFeature) {
        face.clear();

        if let Some(normal) = self.normal() {
            face.set_feature_id(id);

            match id.unwrap_face() {
                0 => {
                    face.push(self.a, FeatureId::Vertex(0));
                    face.push(self.b, FeatureId::Vertex(1));
                    face.push(self.c, FeatureId::Vertex(2));
                    face.push_edge_feature_id(FeatureId::Edge(0));
                    face.push_edge_feature_id(FeatureId::Edge(1));
                    face.push_edge_feature_id(FeatureId::Edge(2));
                    face.set_normal(normal);
                }
                1 => {
                    face.push(self.a, FeatureId::Vertex(0));
                    face.push(self.c, FeatureId::Vertex(2));
                    face.push(self.b, FeatureId::Vertex(1));
                    face.push_edge_feature_id(FeatureId::Edge(2));
                    face.push_edge_feature_id(FeatureId::Edge(1));
                    face.push_edge_feature_id(FeatureId::Edge(0));
                    face.set_normal(-normal);
                }
                _ => unreachable!(),
            }

            face.recompute_edge_normals();
        } else {
            face.push(self.a, FeatureId::Vertex(0));
            face.set_feature_id(FeatureId::Vertex(0));
        }
    }

    fn support_face_toward(
        &self,
        m: &Isometry<Real>,
        dir: &Unit<Vector<Real>>,
        face: &mut ConvexPolygonalFeature,
    ) {
        let normal = self.scaled_normal();

        if normal.dot(&*dir) >= 0.0 {
            ConvexPolyhedron::face(self, FeatureId::Face(0), face);
        } else {
            ConvexPolyhedron::face(self, FeatureId::Face(1), face);
        }
        face.transform_by(m)
    }

    fn support_feature_toward(
        &self,
        transform: &Isometry<Real>,
        dir: &Unit<Vector<Real>>,
        eps: Real,
        out: &mut ConvexPolygonalFeature,
    ) {
        out.clear();
        let tri = self.transformed(transform);
        let feature = tri.support_feature_id_toward(dir, eps);

        match feature {
            FeatureId::Vertex(_) => {
                let v = tri.vertex(feature);
                out.push(v, feature);
                out.set_feature_id(feature);
            }
            FeatureId::Edge(_) => {
                let (a, b, fa, fb) = tri.edge(feature);
                out.push(a, fa);
                out.push(b, fb);
                out.push_edge_feature_id(feature);
                out.set_feature_id(feature);
            }
            FeatureId::Face(_) => tri.face(feature, out),
            _ => unreachable!(),
        }
    }

    fn support_feature_id_toward(&self, local_dir: &Unit<Vector<Real>>) -> FeatureId {
        self.support_feature_id_toward(local_dir, na::convert::<f64, Real>(f64::consts::PI / 180.0))
    }
}
*/

#[cfg(feature = "dim3")]
#[cfg(test)]
mod test {
    use crate::shape::Triangle;
    use na::Point3;

    #[test]
    fn test_triangle_area() {
        let pa = Point3::new(0.0, 5.0, 0.0);
        let pb = Point3::new(0.0, 0.0, 0.0);
        let pc = Point3::new(0.0, 0.0, 4.0);

        assert!(relative_eq!(Triangle::new(pa, pb, pc).area(), 10.0));
    }
}
