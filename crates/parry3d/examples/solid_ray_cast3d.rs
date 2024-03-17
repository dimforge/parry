use parry3d::math::{GlamVectorOps, Isometry, Point, PointOps, Vector};
use parry3d::query::{Ray, RayCast};
use parry3d::shape::Cuboid;

fn main() {
    let cuboid = Cuboid::new(Vector::new(1.0, 2.0, 1.0));
    let ray_inside = Ray::new(Point::origin(), Vector::y());
    let ray_miss = Ray::new(Point::new(2.0, 2.0, 2.0), Vector::new(1.0, 1.0, 1.0));

    // Solid cast.
    assert_eq!(
        cuboid
            .cast_ray(&Isometry::identity(), &ray_inside, std::f32::MAX, true)
            .unwrap(),
        0.0
    );

    // Non-solid cast.
    assert_eq!(
        cuboid
            .cast_ray(&Isometry::identity(), &ray_inside, std::f32::MAX, false)
            .unwrap(),
        2.0
    );

    // The other ray does not intersect this shape.
    assert!(cuboid
        .cast_ray(&Isometry::identity(), &ray_miss, std::f32::MAX, false)
        .is_none());
    assert!(cuboid
        .cast_ray(&Isometry::identity(), &ray_miss, std::f32::MAX, true)
        .is_none());
}
