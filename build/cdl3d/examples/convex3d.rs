extern crate nalgebra as na;

use cdl3d::shape::ConvexPolyhedron;
use na::Point3;

fn main() {
    let points = [
        Point3::new(0.0f32, 0.0, 1.0),
        Point3::new(0.0, 0.0, -1.0),
        Point3::new(0.0, 1.0, 0.0),
        Point3::new(0.0, -1.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(-1.0, 0.0, 0.0),
        Point3::new(0.0, 0.0, 0.0),
    ];

    let convex = ConvexPolyhedron::from_convex_hull(&points).expect("Invalid convex shape.");
    convex.check_geometry();
}
