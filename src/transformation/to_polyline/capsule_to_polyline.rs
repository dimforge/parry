use crate::math::Real;
use crate::shape::Capsule;
use crate::transformation::utils;
use na::{self, Point2, RealField, Vector2};

impl Capsule {
    /// Discretize the boundary of this capsule as a polygonal line.
    pub fn to_polyline(&self, nsubdiv: u32) -> Vec<Point2<Real>> {
        let pi = Real::pi();
        let dtheta = pi / (nsubdiv as Real);

        let mut points: Vec<Point2<Real>> = Vec::with_capacity(nsubdiv as usize);

        utils::push_xy_arc(self.radius, nsubdiv, dtheta, &mut points);

        let npoints = points.len();

        for i in 0..npoints {
            let new_point = points[i] + Vector2::new(na::zero(), self.half_height());

            points.push(-new_point);
            points[i] = new_point;
        }

        utils::transformed(points, self.canonical_transform())
    }
}
