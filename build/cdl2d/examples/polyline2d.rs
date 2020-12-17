extern crate nalgebra as na;

use cdl2d::shape::Polyline;
use na::Point2;

fn main() {
    let points = vec![
        Point2::new(0.0, 1.0),
        Point2::new(-1.0, -1.0),
        Point2::new(0.0, -0.5),
        Point2::new(1.0, -1.0),
    ];

    let indices = vec![
        Point2::new(0, 1),
        Point2::new(1, 2),
        Point2::new(2, 3),
        Point2::new(3, 0), // This forms a loop.
    ];

    // Build the polyline.
    let polyline = Polyline::new(points, Some(indices));

    assert_eq!(polyline.vertices().len(), 4);
}
