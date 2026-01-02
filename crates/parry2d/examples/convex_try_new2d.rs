use parry2d::math::Vector;
use parry2d::shape::ConvexPolygon;

fn main() {
    let points = vec![
        Vector::new(-1.0, 1.0),
        Vector::new(-0.5, -0.5),
        Vector::new(0.5, -0.5),
        Vector::new(1.0, 1.0),
    ];

    let convex = ConvexPolygon::from_convex_polyline(points).expect("Invalid convex polygon.");
    assert!(convex.points().len() == 4);
}
