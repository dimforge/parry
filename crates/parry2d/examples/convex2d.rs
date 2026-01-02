extern crate num_traits as num;

use parry2d::math::Vector;
use parry2d::shape::ConvexPolygon;

fn main() {
    let points = [
        Vector::new(-1.0, 1.0),
        Vector::new(-0.5, -0.5),
        Vector::new(0.0, 0.5),
        Vector::new(0.5, -0.5),
        Vector::new(1.0, 1.0),
    ];

    let convex = ConvexPolygon::from_convex_hull(&points).expect("Invalid convex polygon.");
    assert!(convex.points().len() == 4);
}
