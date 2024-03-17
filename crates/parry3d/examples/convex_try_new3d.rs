use parry3d::math::Point;
use parry3d::shape::ConvexPolyhedron;

fn main() {
    let points = vec![
        Point::new(0.0f32, 0.0, 1.0),
        Point::new(0.0, 0.0, -1.0),
        Point::new(0.0, 1.0, 0.0),
        Point::new(0.0, -1.0, 0.0),
        Point::new(1.0, 0.0, 0.0),
        Point::new(-1.0, 0.0, 0.0),
    ];

    let indices = vec![
        [0, 4, 2],
        [0, 3, 4],
        [5, 0, 2],
        [5, 3, 0],
        [1, 5, 2],
        [1, 3, 5],
        [4, 1, 2],
        [4, 3, 1],
    ];

    let convex =
        ConvexPolyhedron::from_convex_mesh(points, &indices).expect("Invalid convex shape.");
    convex.check_geometry();
}
