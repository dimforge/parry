extern crate nalgebra as na;
extern crate num_traits as num;

use na::Point2;
use parry2d::shape::ConvexPolygon;

fn main() {
    let points = [
        Point2::new(-1.0, 1.0),
        Point2::new(-0.5, -0.5),
        Point2::new(0.0, 0.5),
        Point2::new(0.5, -0.5),
        Point2::new(1.0, 1.0),
    ];

    let convex = ConvexPolygon::from_convex_hull(&points).expect("Invalid convex polygon.");
    assert!(convex.points().len() == 4);
}
