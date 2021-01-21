use crate::math::{Point, Real};
use crate::query::gjk::{self, CSOPoint};
use crate::query::{PointQuery, PointQueryWithLocation};
use crate::shape::{Segment, SegmentPointLocation, Triangle, TrianglePointLocation};

/// A simplex of dimension up to 2 using Voronoï regions for computing point projections.
#[derive(Clone, Debug)]
pub struct VoronoiSimplex {
    prev_vertices: [usize; 3],
    prev_dim: usize,
    prev_proj: [Real; 2],

    vertices: [CSOPoint; 3],
    proj: [Real; 2],
    dim: usize,
}

impl VoronoiSimplex {
    /// Crates a new empty simplex.
    pub fn new() -> VoronoiSimplex {
        VoronoiSimplex {
            prev_vertices: [0, 1, 2],
            prev_proj: [0.0; 2],
            prev_dim: 0,
            vertices: [CSOPoint::origin(); 3],
            proj: [0.0; 2],
            dim: 0,
        }
    }

    /// Swap two vertices of this simplex.
    pub fn swap(&mut self, i1: usize, i2: usize) {
        self.vertices.swap(i1, i2);
        self.prev_vertices.swap(i1, i2);
    }

    /// Resets this simplex to a single point.
    pub fn reset(&mut self, pt: CSOPoint) {
        self.prev_dim = 0;
        self.dim = 0;
        self.vertices[0] = pt;
    }

    /// Add a point to this simplex.
    pub fn add_point(&mut self, pt: CSOPoint) -> bool {
        self.prev_dim = self.dim;
        self.prev_proj = self.proj;
        self.prev_vertices = [0, 1, 2];

        for i in 0..self.dim + 1 {
            if (self.vertices[i].point - pt.point).norm_squared() < gjk::eps_tol() {
                return false;
            }
        }

        self.dim += 1;
        self.vertices[self.dim] = pt;
        return true;
    }

    /// Retrieves the barycentric coordinate associated to the `i`-th by the last call to `project_origin_and_reduce`.
    pub fn proj_coord(&self, i: usize) -> Real {
        assert!(i <= self.dim, "Index out of bounds.");
        self.proj[i]
    }

    /// The i-th point of this simplex.
    pub fn point(&self, i: usize) -> &CSOPoint {
        assert!(i <= self.dim, "Index out of bounds.");
        &self.vertices[i]
    }

    /// Retrieves the barycentric coordinate associated to the `i`-th before the last call to `project_origin_and_reduce`.
    pub fn prev_proj_coord(&self, i: usize) -> Real {
        assert!(i <= self.prev_dim, "Index out of bounds.");
        self.prev_proj[i]
    }

    /// The i-th point of the simplex before the last call to `projet_origin_and_reduce`.
    pub fn prev_point(&self, i: usize) -> &CSOPoint {
        assert!(i <= self.prev_dim, "Index out of bounds.");
        &self.vertices[self.prev_vertices[i]]
    }

    /// Projets the origin on the boundary of this simplex and reduces `self` the smallest subsimplex containing the origin.
    ///
    /// Retruns the result of the projection or Point::origin() if the origin lies inside of the simplex.
    /// The state of the samplex before projection is saved, and can be retrieved using the methods prefixed
    /// by `prev_`.
    pub fn project_origin_and_reduce(&mut self) -> Point<Real> {
        if self.dim == 0 {
            self.proj[0] = 1.0;
            self.vertices[0].point
        } else if self.dim == 1 {
            // FIXME: NLL
            let (proj, location) = {
                let seg = Segment::new(self.vertices[0].point, self.vertices[1].point);
                seg.project_local_point_and_get_location(&Point::<Real>::origin(), true)
            };

            match location {
                SegmentPointLocation::OnVertex(0) => {
                    self.proj[0] = 1.0;
                    self.dim = 0;
                }
                SegmentPointLocation::OnVertex(1) => {
                    self.proj[0] = 1.0;
                    self.swap(0, 1);
                    self.dim = 0;
                }
                SegmentPointLocation::OnEdge(coords) => {
                    self.proj = coords;
                }
                _ => unreachable!(),
            }

            proj.point
        } else {
            assert!(self.dim == 2);
            // FIXME: NLL
            let (proj, location) = {
                let tri = Triangle::new(
                    self.vertices[0].point,
                    self.vertices[1].point,
                    self.vertices[2].point,
                );
                tri.project_local_point_and_get_location(&Point::<Real>::origin(), true)
            };

            match location {
                TrianglePointLocation::OnVertex(i) => {
                    self.swap(0, i as usize);
                    self.proj[0] = 1.0;
                    self.dim = 0;
                }
                TrianglePointLocation::OnEdge(0, coords) => {
                    self.proj = coords;
                    self.dim = 1;
                }
                TrianglePointLocation::OnEdge(1, coords) => {
                    self.swap(0, 2);
                    self.proj[0] = coords[1];
                    self.proj[1] = coords[0];
                    self.dim = 1;
                }
                TrianglePointLocation::OnEdge(2, coords) => {
                    self.swap(1, 2);
                    self.proj = coords;
                    self.dim = 1;
                }
                _ => {}
            }

            proj.point
        }
    }

    /// Compute the projection of the origin on the boundary of this simplex.
    pub fn project_origin(&mut self) -> Point<Real> {
        if self.dim == 0 {
            self.vertices[0].point
        } else if self.dim == 1 {
            let seg = Segment::new(self.vertices[0].point, self.vertices[1].point);
            seg.project_local_point(&Point::<Real>::origin(), true)
                .point
        } else {
            assert!(self.dim == 2);
            let tri = Triangle::new(
                self.vertices[0].point,
                self.vertices[1].point,
                self.vertices[2].point,
            );
            tri.project_local_point(&Point::<Real>::origin(), true)
                .point
        }
    }

    /// Tests if the given point is already a vertex of this simplex.
    pub fn contains_point(&self, pt: &Point<Real>) -> bool {
        for i in 0..self.dim + 1 {
            if self.vertices[i].point == *pt {
                return true;
            }
        }

        false
    }

    /// The dimension of the smallest subspace that can contain this simplex.
    pub fn dimension(&self) -> usize {
        self.dim
    }

    /// The dimension of the simplex before the last call to `project_origin_and_reduce`.
    pub fn prev_dimension(&self) -> usize {
        self.prev_dim
    }

    /// The maximum squared length of the vertices of this simplex.
    pub fn max_sq_len(&self) -> Real {
        let mut max_sq_len = 0.0;

        for i in 0..self.dim + 1 {
            let norm = self.vertices[i].point.coords.norm_squared();

            if norm > max_sq_len {
                max_sq_len = norm
            }
        }

        max_sq_len
    }

    /// Apply a function to all the vertices of this simplex.
    pub fn modify_pnts(&mut self, f: &dyn Fn(&mut CSOPoint)) {
        for i in 0..self.dim + 1 {
            f(&mut self.vertices[i])
        }
    }
}
