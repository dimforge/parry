extern crate nalgebra as na;

use parry3d::math::Point;
use parry3d::shape::Polyline;

fn main() {
    let points = vec![
        Point::new(0.0, 1.0, 0.0),
        Point::new(-1.0, -1.0, 1.0),
        Point::new(0.0, -0.5, 0.0),
        Point::new(1.0, -1.0, -1.0),
        Point::new(0.0, 1.0, 0.0), // This forms a loop.
    ];

    // Build the polyline.
    let polyline = Polyline::new(points, None);

    assert!(polyline.vertices().len() == 5);
}
