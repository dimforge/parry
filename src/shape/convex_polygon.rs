use crate::math::{Point, Real, Vector};
use crate::shape::{FeatureId, PackedFeatureId, PolygonalFeature, PolygonalFeatureMap, SupportMap};
use crate::utils;
use na::{self, ComplexField, RealField, Unit};

/// A 2D convex polygon.
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Deserialize, rkyv::Serialize),
    archive(check_bytes)
)]
#[derive(Clone, Debug)]
pub struct ConvexPolygon {
    points: Vec<Point<Real>>,
    normals: Vec<Unit<Vector<Real>>>,
}

impl ConvexPolygon {
    /// Creates a new 2D convex polygon from an arbitrary set of points.
    ///
    /// This explicitly computes the convex hull of the given set of points.
    /// Returns `None` if the convex hull computation failed.
    pub fn from_convex_hull(points: &[Point<Real>]) -> Option<Self> {
        let vertices = crate::transformation::convex_hull(points);
        Self::from_convex_polyline(vertices)
    }

    /// Creates a new 2D convex polygon from a set of points assumed to
    /// describe a counter-clockwise convex polyline.
    ///
    /// Convexity of the input polyline is not checked.
    /// Returns `None` if all points form an almost flat line.
    pub fn from_convex_polyline(mut points: Vec<Point<Real>>) -> Option<Self> {
        let eps = ComplexField::sqrt(crate::math::DEFAULT_EPSILON);
        let mut normals = Vec::with_capacity(points.len());

        // First, compute all normals.
        for i1 in 0..points.len() {
            let i2 = (i1 + 1) % points.len();
            normals.push(utils::ccw_face_normal([&points[i1], &points[i2]])?);
        }

        let mut nremoved = 0;
        // See if the first vertex must be removed.
        if normals[0].dot(&*normals[normals.len() - 1]) > 1.0 - eps {
            nremoved = 1;
        }

        // Second, find vertices that can be removed because
        // of collinearity of adjascent faces.
        for i2 in 1..points.len() {
            let i1 = i2 - 1;
            if normals[i1].dot(&*normals[i2]) > 1.0 - eps {
                // Remove
                nremoved += 1;
            } else {
                points[i2 - nremoved] = points[i2];
                normals[i2 - nremoved] = normals[i2];
            }
        }

        let new_length = points.len() - nremoved;
        points.truncate(new_length);
        normals.truncate(new_length);

        if !points.is_empty() {
            Some(ConvexPolygon { points, normals })
        } else {
            None
        }
    }

    /// The vertices of this convex polygon.
    #[inline]
    pub fn points(&self) -> &[Point<Real>] {
        &self.points
    }

    /// The normals of the edges of this convex polygon.
    #[inline]
    pub fn normals(&self) -> &[Unit<Vector<Real>>] {
        &self.normals
    }

    /// Computes a scaled version of this convex polygon.
    ///
    /// Returns `None` if the result had degenerate normals (for example if
    /// the scaling factor along one axis is zero).
    pub fn scaled(mut self, scale: &Vector<Real>) -> Option<Self> {
        self.points
            .iter_mut()
            .for_each(|pt| pt.coords.component_mul_assign(scale));

        for n in &mut self.normals {
            *n = Unit::try_new(n.component_mul(scale), 0.0)?;
        }

        Some(self)
    }

    /// Returns a mitered offset of the polygon.
    ///
    /// # Arguments
    ///
    /// * `amount` - size of the inflation. Each edge is moved outwards by this
    ///   amount.
    ///
    /// # Panics
    ///
    /// Panics if `amount` is not a non-negative finite number.
    pub fn offsetted(&self, amount: Real) -> Self {
        if !amount.is_finite() || amount < 0. {
            panic!(
                "Offset amount must be a non-negative finite number, got {}.",
                amount
            );
        }

        let mut points = Vec::with_capacity(self.points.len());
        let normals = self.normals.clone();

        for i2 in 0..self.points.len() {
            let i1 = if i2 == 0 {
                self.points.len() - 1
            } else {
                i2 - 1
            };
            let normal_a = normals[i1];
            let direction = normal_a.into_inner() + normals[i2].into_inner();
            points.push(self.points[i2] + (amount / direction.dot(&normal_a)) * direction);
        }

        ConvexPolygon { points, normals }
    }

    /// Get the ID of the feature with a normal that maximizes the dot product with `local_dir`.
    pub fn support_feature_id_toward(&self, local_dir: &Unit<Vector<Real>>) -> FeatureId {
        let eps: Real = Real::pi() / 180.0;
        let ceps = ComplexField::cos(eps);

        // Check faces.
        for i in 0..self.normals.len() {
            let normal = &self.normals[i];

            if normal.dot(local_dir.as_ref()) >= ceps {
                return FeatureId::Face(i as u32);
            }
        }

        // Support vertex.
        FeatureId::Vertex(
            utils::point_cloud_support_point_id(local_dir.as_ref(), &self.points) as u32,
        )
    }

    /// The normal of the given feature.
    pub fn feature_normal(&self, feature: FeatureId) -> Option<Unit<Vector<Real>>> {
        match feature {
            FeatureId::Face(id) => Some(self.normals[id as usize]),
            FeatureId::Vertex(id2) => {
                let id1 = if id2 == 0 {
                    self.normals.len() - 1
                } else {
                    id2 as usize - 1
                };
                Some(Unit::new_normalize(
                    *self.normals[id1] + *self.normals[id2 as usize],
                ))
            }
            _ => None,
        }
    }
}

impl SupportMap for ConvexPolygon {
    #[inline]
    fn local_support_point(&self, dir: &Vector<Real>) -> Point<Real> {
        utils::point_cloud_support_point(dir, self.points())
    }
}

impl PolygonalFeatureMap for ConvexPolygon {
    fn local_support_feature(&self, dir: &Unit<Vector<Real>>, out_feature: &mut PolygonalFeature) {
        let cuboid = crate::shape::Cuboid::new(self.points[2].coords);
        cuboid.local_support_feature(dir, out_feature);
        let mut best_face = 0;
        let mut max_dot = self.normals[0].dot(dir);

        for i in 1..self.normals.len() {
            let dot = self.normals[i].dot(dir);

            if dot > max_dot {
                max_dot = dot;
                best_face = i;
            }
        }

        let i1 = best_face;
        let i2 = (best_face + 1) % self.points.len();
        *out_feature = PolygonalFeature {
            vertices: [self.points[i1], self.points[i2]],
            vids: PackedFeatureId::vertices([i1 as u32 * 2, i2 as u32 * 2]),
            fid: PackedFeatureId::face(i1 as u32 * 2 + 1),
            num_vertices: 2,
        };
    }
}

/*
impl ConvexPolyhedron for ConvexPolygon {
    fn vertex(&self, id: FeatureId) -> Point<Real> {
        self.points[id.unwrap_vertex() as usize]
    }

    fn face(&self, id: FeatureId, out: &mut ConvexPolygonalFeature) {
        out.clear();

        let ia = id.unwrap_face() as usize;
        let ib = (ia + 1) % self.points.len();
        out.push(self.points[ia], FeatureId::Vertex(ia as u32));
        out.push(self.points[ib], FeatureId::Vertex(ib as u32));

        out.set_normal(self.normals[ia as usize]);
        out.set_feature_id(FeatureId::Face(ia as u32));
    }



    fn support_face_toward(
        &self,
        m: &Isometry<Real>,
        dir: &Unit<Vector<Real>>,
        out: &mut ConvexPolygonalFeature,
    ) {
        let ls_dir = m.inverse_transform_vector(dir);
        let mut best_face = 0;
        let mut max_dot = self.normals[0].dot(&ls_dir);

        for i in 1..self.points.len() {
            let dot = self.normals[i].dot(&ls_dir);

            if dot > max_dot {
                max_dot = dot;
                best_face = i;
            }
        }

        self.face(FeatureId::Face(best_face as u32), out);
        out.transform_by(m);
    }

    fn support_feature_toward(
        &self,
        transform: &Isometry<Real>,
        dir: &Unit<Vector<Real>>,
        _angle: Real,
        out: &mut ConvexPolygonalFeature,
    ) {
        out.clear();
        // FIXME: actualy find the support feature.
        self.support_face_toward(transform, dir, out)
    }

    fn support_feature_id_toward(&self, local_dir: &Unit<Vector<Real>>) -> FeatureId {
        let eps: Real = na::convert::<f64, Real>(f64::consts::PI / 180.0);
        let ceps = ComplexField::cos(eps);

        // Check faces.
        for i in 0..self.normals.len() {
            let normal = &self.normals[i];

            if normal.dot(local_dir.as_ref()) >= ceps {
                return FeatureId::Face(i as u32);
            }
        }

        // Support vertex.
        FeatureId::Vertex(
            utils::point_cloud_support_point_id(local_dir.as_ref(), &self.points) as u32,
        )
    }
}
*/

#[cfg(feature = "dim2")]
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_dilation() {
        let polygon = ConvexPolygon::from_convex_polyline(vec![
            Point::new(1., 0.),
            Point::new(-1., 0.),
            Point::new(0., -1.),
        ])
        .unwrap();

        let offsetted = polygon.offsetted(0.5);
        let expected = vec![
            Point::new(2.207, 0.5),
            Point::new(-2.207, 0.5),
            Point::new(0., -1.707),
        ];

        assert_eq!(offsetted.points().len(), 3);
        assert!(offsetted
            .points()
            .iter()
            .zip(expected.iter())
            .all(|(a, b)| (a.coords - b.coords).magnitude() < 0.001));
    }
}
