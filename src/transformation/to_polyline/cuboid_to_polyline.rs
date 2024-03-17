use crate::math::*;
use crate::shape::Cuboid;
use crate::transformation::utils;
use na::{self};

impl Cuboid {
    /// Discretize the boundary of this cuboid as a polygonal line.
    pub fn to_polyline(&self) -> Vec<Point> {
        utils::scaled(unit_rectangle(), self.half_extents * 2.0)
    }
}

/// The contour of a unit cuboid lying on the x-y plane.
fn unit_rectangle() -> Vec<Point> {
    vec![
        Point::new(-0.5, -0.5),
        Point::new(0.5, -0.5),
        Point::new(0.5, 0.5),
        Point::new(-0.5, 0.5),
    ]
}
