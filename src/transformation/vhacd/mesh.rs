use crate::shape::Tetrahedron;
use bitflags::_core::cmp::min;
use na::{Point3, Vector3};
use std::convert::TryFrom;

#[derive(Copy, Clone, Debug)]
pub enum Axis {
    X = 0,
    Y,
    Z,
}

impl Axis {
    pub fn try_from(val: usize) -> Option<Self> {
        match val {
            0 => Some(Axis::X),
            1 => Some(Axis::Y),
            2 => Some(Axis::Z),
            _ => None,
        }
    }
}

#[derive(Copy, Clone, Debug)]
pub struct Plane {
    pub abc: Vector3<Real>,
    pub d: Real,
    pub axis: Axis,
    pub index: u32,
}

pub struct Mesh {
    pub points: Vec<Point3<Real>>,
    pub triangles: Vec<Point3<u32>>,
    pub min_bb: Point3<Real>,
    pub max_bb: Point3<Real>,
    pub diag: Real,
    pub center: Point3<Real>,
}

impl Mesh {
    pub fn new() -> Self {
        Self {
            points: Vec::new(),
            triangles: Vec::new(),
            min_bb: Point3::origin(),
            max_bb: Point3::origin(),
            diag: 0.0,
            center: Point3::origin(),
        }
    }

    pub(crate) fn compute_convex_hull(&mut self, points: &[Point3<Real>]) {
        // Set `self` to the convex hull of `points`.
        unimplemented!()
    }

    // Update the mesh center and its bounding-box.
    pub(crate) fn compute_center(&mut self) {
        if !self.points.is_empty() {
            // FIXME: use the ncollide center function instead.
            self.center = crate::utils::center(&self.points);

            let aabb = crate::bounding_volume::local_point_cloud_aabb(&self.points);
            self.min_bb = aabb.mins;
            self.max_bb = aabb.maxs;
        }
    }

    /// Computes the volume of this mesh, assuming it is convex.
    // TODO: rename this "compute_convex_volume".
    pub(crate) fn compute_volume(&self) -> Real {
        let num_points = self.points.len();
        let num_triangles = self.triangles.len();

        if num_points == 0 || num_triangles == 0 {
            return 0.0;
        }

        let barycenter = crate::utils::center(&self.points);
        let mut total_volume = 0.0;

        for tri in &self.triangles {
            let a = self.points[tri.x as usize];
            let b = self.points[tri.y as usize];
            let c = self.points[tri.z as usize];
            total_volume += Tetrahedron::new(a, b, c, barycenter).volume();
        }

        total_volume
    }

    /// Split the vertices of this mesh into two parts, depending on what side of the given plan they lie.
    ///
    /// Points located exactly on the plane boundary are pushed to both positive and negative parts.
    pub(crate) fn clip(
        &self,
        plane: &Plane,
        positive_part: &mut Vec<Point3<Real>>,
        negative_part: &mut Vec<Point3<Real>>,
    ) {
        for pt in &self.points {
            let d = plane.abc.dot(&pt.coords) + plane.d;

            if d > 0.0 {
                positive_part.push(*pt);
            } else if d < 0.0 {
                negative_part.push(*pt);
            } else {
                positive_part.push(*pt);
                negative_part.push(*pt);
            }
        }
    }

    /// Tests if the given point is inside of this mesh, assuming this mesh is convex.
    pub(super) fn is_inside(&self, pt: &Point3<Real>) -> bool {
        if self.points.is_empty() || self.triangles.is_empty() {
            return false;
        }

        for tri in &self.triangles {
            let ver0 = self.points[tri[0] as usize];
            let ver1 = self.points[tri[1] as usize];
            let ver2 = self.points[tri[2] as usize];
            let volume = Tetrahedron::new(ver0, ver1, ver2, *pt).signed_volume();

            if volume < 0.0 {
                return false;
            }
        }

        true
    }

    /// Updates the diagonal of the bounding box of this mesh.
    fn compute_diag_bb(&mut self) {
        if self.points.is_empty() {
            return;
        }

        let aabb = crate::bounding_volume::local_point_cloud_aabb(&self.points);
        self.diag = aabb.extents().norm();
    }
}

/*
void Mesh::ComputeConvexHull(const double* const pts, const size_t nPts)
{
ResizePoints(0);
ResizeTriangles(0);
btConvexHullComputer ch;
ch.compute(pts, 3 * sizeof(double), (int32_t)nPts, -1.0, -1.0);
for (int32_t v = 0; v < ch.vertices.size(); v++)
{
AddPoint(Vec3<double>(ch.vertices[v].getX(), ch.vertices[v].getY(), ch.vertices[v].getZ()));
}
const int32_t nt = ch.faces.size();
for (int32_t t = 0; t < nt; ++t)
{
const btConvexHullComputer::Edge* sourceEdge = &(ch.edges[ch.faces[t]]);
int32_t a = sourceEdge->getSourceVertex();
int32_t b = sourceEdge->getTargetVertex();
const btConvexHullComputer::Edge* edge = sourceEdge->getNextEdgeOfFace();
int32_t c = edge->getTargetVertex();
while (c != a)
{
AddTriangle(Vec3<int32_t>(a, b, c));
edge = edge->getNextEdgeOfFace();
b = c;
c = edge->getTargetVertex();
}
}
}
*/
