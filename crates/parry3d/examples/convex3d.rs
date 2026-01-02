use parry3d::math::Vector;
use parry3d::shape::ConvexPolyhedron;

fn main() {
    let points = [
        Vector::new(0.0f32, 0.0, 1.0),
        Vector::new(0.0, 0.0, -1.0),
        Vector::new(0.0, 1.0, 0.0),
        Vector::new(0.0, -1.0, 0.0),
        Vector::new(1.0, 0.0, 0.0),
        Vector::new(-1.0, 0.0, 0.0),
        Vector::ZERO,
    ];

    let convex = ConvexPolyhedron::from_convex_hull(&points).expect("Invalid convex shape.");
    convex.check_geometry();
}
