use na::{self, Isometry2, Point2, Vector2};
use parry2d::math::Real;
use parry2d::query::{Ray, RayCast};
use parry2d::shape::{ConvexPolygon, Segment, Triangle};

#[test]
fn issue_178_parallel_raycast() {
    let m1 = Isometry2::identity();
    let ray = Ray::new(Point2::new(0.0, 0.0), Vector2::new(0.0, 1.0));
    let seg = Segment::new(Point2::new(2.0, 1.0), Point2::new(2.0, 0.0));

    let cast = seg.cast_ray(&m1, &ray, std::f32::MAX, true);
    assert!(cast.is_none());
}

#[test]
fn parallel_raycast() {
    let m1 = Isometry2::identity();
    let ray = Ray::new(Point2::new(0.0, 0.0), Vector2::new(0.0, 1.0));
    let seg = Segment::new(Point2::new(2.0, 1.0), Point2::new(2.0, -1.0));

    let cast = seg.cast_ray(&m1, &ray, std::f32::MAX, true);
    assert!(cast.is_none());
}

#[test]
fn collinear_raycast_starting_on_segment() {
    let m1 = Isometry2::identity();
    let ray = Ray::new(Point2::new(0.0, 0.0), Vector2::new(0.0, 1.0));
    let seg = Segment::new(Point2::new(0.0, 1.0), Point2::new(0.0, -1.0));

    let cast = seg.cast_ray(&m1, &ray, std::f32::MAX, true);
    assert_eq!(cast, Some(0.0));
}

#[test]
fn collinear_raycast_starting_below_segment() {
    let m1 = Isometry2::identity();
    let ray = Ray::new(Point2::new(0.0, -2.0), Vector2::new(0.0, 1.0));
    let seg = Segment::new(Point2::new(0.0, 1.0), Point2::new(0.0, -1.0));

    let cast = seg.cast_ray(&m1, &ray, std::f32::MAX, true);
    assert_eq!(cast, Some(1.0));
}

#[test]
fn collinear_raycast_starting_above_segment() {
    let m1 = Isometry2::identity();
    let ray = Ray::new(Point2::new(0.0, 2.0), Vector2::new(0.0, 1.0));
    let seg = Segment::new(Point2::new(0.0, 1.0), Point2::new(0.0, -1.0));

    let cast = seg.cast_ray(&m1, &ray, std::f32::MAX, true);
    assert_eq!(cast, None);
}

#[test]
fn perpendicular_raycast_starting_behind_segment() {
    let segment = Segment::new(Point2::new(0.0f32, -10.0), Point2::new(0.0, 10.0));
    let ray = Ray::new(Point2::new(-1.0, 0.0), Vector2::new(1.0, 0.0));
    assert!(segment.intersects_local_ray(&ray, std::f32::MAX));
}

#[test]
fn perpendicular_raycast_starting_in_front_of_segment() {
    let segment = Segment::new(Point2::new(0.0f32, -10.0), Point2::new(0.0, 10.0));
    let ray = Ray::new(Point2::new(1.0, 0.0), Vector2::new(1.0, 0.0));
    assert!(!segment.intersects_local_ray(&ray, std::f32::MAX));
}

#[test]
fn perpendicular_raycast_starting_on_segment() {
    let segment = Segment::new(Point2::new(0.0f32, -10.0), Point2::new(0.0, 10.0));
    let ray = Ray::new(Point2::new(0.0, 3.0), Vector2::new(1.0, 0.0));

    let cast = segment.cast_local_ray(&ray, std::f32::MAX, true);
    assert_eq!(cast, Some(0.0));
}

#[test]
fn perpendicular_raycast_starting_above_segment() {
    let segment = Segment::new(Point2::new(0.0f32, -10.0), Point2::new(0.0, 10.0));
    let ray = Ray::new(Point2::new(0.0, 11.0), Vector2::new(1.0, 0.0));
    assert!(!segment.intersects_local_ray(&ray, std::f32::MAX));
}

#[test]
fn perpendicular_raycast_starting_below_segment() {
    let segment = Segment::new(Point2::new(0.0f32, -10.0), Point2::new(0.0, 10.0));
    let ray = Ray::new(Point2::new(0.0, -11.0), Vector2::new(1.0, 0.0));
    assert!(!segment.intersects_local_ray(&ray, std::f32::MAX));
}

#[test]
fn raycast_starting_outside_of_triangle() {
    let triangle = Triangle::new(
        Point2::new(0.0f32, -10.0),
        Point2::new(0.0, 10.0),
        Point2::new(10.0, 0.0),
    );
    let ray = Ray::new(Point2::new(-10.0, 0.0), Vector2::new(1.0, 0.0));
    let intersect = triangle
        .cast_local_ray_and_get_normal(&ray, std::f32::MAX, true)
        .expect("No intersection");

    assert_ne!(intersect.time_of_impact, 0.0);
}

#[test]
fn raycast_starting_inside_of_triangle() {
    let triangle = Triangle::new(
        Point2::new(0.0f32, -10.0),
        Point2::new(0.0, 10.0),
        Point2::new(10.0, 0.0),
    );
    let ray = Ray::new(Point2::new(2.0, 0.0), Vector2::new(1.0, 0.0));
    let intersect = triangle
        .cast_local_ray_and_get_normal(&ray, std::f32::MAX, true)
        .expect("No intersection");

    assert_eq!(intersect.time_of_impact, 0.0);
}

#[test]
fn raycast_starting_on_edge_of_triangle() {
    let triangle = Triangle::new(
        Point2::new(0.0f32, -10.0),
        Point2::new(0.0, 10.0),
        Point2::new(10.0, 0.0),
    );
    let ray = Ray::new(Point2::new(0.0, 0.0), Vector2::new(1.0, 0.0));
    let intersect = triangle
        .cast_local_ray_and_get_normal(&ray, std::f32::MAX, true)
        .expect("No intersection");

    assert_eq!(intersect.time_of_impact, 0.0);
}

///    Ray Target
///    +
/// 3  |     Ray moved up each iteration
/// |  v     ^
/// 2  *--+  |
/// |  |  |  +
/// 1  +--+  *<----+Ray origin
/// |
/// 0--1--2--3
///
/// Tests the accuracy of raycaster collision detection against a `ConvexPolygon`.
#[test]
fn convexpoly_raycast_fuzz() {
    let vertices: Vec<Point2<Real>> = vec![[2, 1], [2, 2], [1, 2], [1, 1]]
        .into_iter()
        .map(|[x, y]| Point2::new(x as Real, y as Real))
        .collect();

    let square = ConvexPolygon::from_convex_polyline(vertices).unwrap();

    let test_raycast = |ray_origin: Point2<Real>, ray_look_at: Point2<Real>| -> Option<Real> {
        let ray_angle = ray_look_at - ray_origin;

        square.cast_ray(
            &Isometry2::identity(),
            &Ray::new(ray_origin, ray_angle.normalize()),
            Real::MAX,
            true,
        )
    };

    for i in 0..10_000 {
        let ray_origin = Point2::new(3., 1. + (i as Real * 0.0001));
        let ray_look_at = Point2::new(0., 2.);
        let collision = test_raycast(ray_origin, ray_look_at);
        let eps = 1.0e-5;

        match collision {
            Some(distance) if distance >= 1.0 - eps && distance < (2.0 as Real).sqrt() => (),
            Some(distance) if distance >= 2.0 => panic!(
                "Collided with back face instead of front face. Distance: {}",
                distance
            ),
            Some(distance) => panic!("Invalid collision distance: {}", distance),
            None => panic!("Failed to collide with any face"),
        }
    }
}
