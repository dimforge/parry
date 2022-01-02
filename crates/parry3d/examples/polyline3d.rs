extern crate nalgebra as na;

use na::Point3;
use parry3d::shape::Polyline;

fn main() {
    let points = vec![
        Point3::new(0.0, 1.0, 0.0),
        Point3::new(-1.0, -1.0, 1.0),
        Point3::new(0.0, -0.5, 0.0),
        Point3::new(1.0, -1.0, -1.0),
        Point3::new(0.0, 1.0, 0.0), // This forms a loop.
    ];

    // Build the polyline.
    let polyline = Polyline::new(points, None);

    assert!(polyline.vertices().len() == 5);
}
