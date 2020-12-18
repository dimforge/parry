use crate::math::{Isometry, Point, Real, Vector};
use crate::shape::{FeatureId, SupportMap};
use crate::utils;
use na::{self, Unit};
use std::f64;

/// A 2D convex polygon.
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[derive(Clone, Debug)]
pub struct ConvexPolygon {
    points: Vec<Point<Real>>,
    normals: Vec<Unit<Vector<Real>>>,
}

impl ConvexPolygon {
    /// Creates a new 2D convex polygon from an arbitrary set of points.
    ///
    /// This explicitly computes the convex hull of the given set of points. Use
    /// Returns `None` if the convex hull computation failed.
    pub fn try_from_points(points: &[Point<Real>]) -> Option<Self> {
        let mut vertices = crate::transformation::convex_hull(points);
        vertices.reverse(); // FIXME: it is unfortunate to have to do this reverse.

        Self::try_new(vertices)
    }

    /// Creates a new 2D convex polygon from a set of points assumed to describe a counter-clockwise convex polyline.
    ///
    /// Convexity of the input polyline is not checked.
    /// Returns `None` if some consecutive points are identical (or too close to being so).
    pub fn try_new(mut points: Vec<Point<Real>>) -> Option<Self> {
        let eps = crate::math::DEFAULT_EPSILON.sqrt();
        let mut normals = Vec::with_capacity(points.len());

        // First, compute all normals.
        for i1 in 0..points.len() {
            let i2 = (i1 + 1) % points.len();
            normals.push(utils::ccw_face_normal([&points[i1], &points[i2]])?);
        }

        let mut nremoved = 0;
        // See if the first vexrtex must be removed.
        if normals[0].dot(&*normals[normals.len() - 1]) > na::one::<Real>() - eps {
            nremoved = 1;
        }

        // Second, find vertices that can be removed because
        // of collinearity of adjascent faces.
        for i2 in 1..points.len() {
            let i1 = i2 - 1;
            if normals[i1].dot(&*normals[i2]) > na::one::<Real>() - eps {
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

        if points.len() != 0 {
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

    pub fn support_feature_id_toward(&self, local_dir: &Unit<Vector<Real>>) -> FeatureId {
        let eps: Real = na::convert::<f64, Real>(f64::consts::PI / 180.0);
        let ceps = eps.cos();

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

impl SupportMap for ConvexPolygon {
    #[inline]
    fn local_support_point(&self, dir: &Vector<Real>) -> Point<Real> {
        utils::point_cloud_support_point(dir, self.points())
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

    fn feature_normal(&self, feature: FeatureId) -> Unit<Vector<Real>> {
        match feature {
            FeatureId::Face(id) => self.normals[id as usize],
            FeatureId::Vertex(id2) => {
                let id1 = if id2 == 0 {
                    self.normals.len() - 1
                } else {
                    id2 as usize - 1
                };
                Unit::new_normalize(*self.normals[id1] + *self.normals[id2 as usize])
            }
            _ => panic!("Invalid feature ID: {:?}", feature),
        }
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
        let ceps = eps.cos();

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