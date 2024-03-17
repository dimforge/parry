use parry3d::math::Point;
use parry3d::shape::ConvexPolyhedron;

fn main() {
    let points = [
        Point::new(0.0f32, 0.0, 1.0),
        Point::new(0.0, 0.0, -1.0),
        Point::new(0.0, 1.0, 0.0),
        Point::new(0.0, -1.0, 0.0),
        Point::new(1.0, 0.0, 0.0),
        Point::new(-1.0, 0.0, 0.0),
        Point::new(0.0, 0.0, 0.0),
    ];

    let convex = ConvexPolyhedron::from_convex_hull(&points).expect("Invalid convex shape.");
    convex.check_geometry();
}
