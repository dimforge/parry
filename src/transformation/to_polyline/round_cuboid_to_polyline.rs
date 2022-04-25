use crate::math::Real;
use crate::shape::Cuboid;
use crate::transformation::utils;
use na::{self, Point2};

impl Cuboid {
    /// Discretize the boundary of this cuboid as a polygonal line.
    pub fn to_polyline(&self) -> Vec<Point2<Real>> {
        utils::scaled(unit_rectangle(), self.half_extents * 2.0)
    }
}

/// The contour of a unit cuboid lying on the x-y plane.
fn unit_rectangle() -> Vec<Point2<Real>> {
    vec![
        Point2::new(-0.5, -0.5),
        Point2::new(0.5, -0.5),
        Point2::new(0.5, 0.5),
        Point2::new(-0.5, 0.5),
    ]
}
