use parry3d::math::Vector;
use parry3d::shape::Polyline;

fn main() {
    let points = vec![
        Vector::new(0.0, 1.0, 0.0),
        Vector::new(-1.0, -1.0, 1.0),
        Vector::new(0.0, -0.5, 0.0),
        Vector::new(1.0, -1.0, -1.0),
        Vector::new(0.0, 1.0, 0.0), // This forms a loop.
    ];

    // Build the polyline.
    let polyline = Polyline::new(points, None);

    assert!(polyline.vertices().len() == 5);
}
