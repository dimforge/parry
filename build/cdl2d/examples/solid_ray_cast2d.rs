extern crate nalgebra as na;

use cdl2d::query::{Ray, RayCast};
use cdl2d::shape::Cuboid;
use na::{Isometry2, Point2, Vector2};

fn main() {
    let cuboid = Cuboid::new(Vector2::new(1.0, 2.0));
    let ray_inside = Ray::new(Point2::origin(), Vector2::y());
    let ray_miss = Ray::new(Point2::new(2.0, 2.0), Vector2::new(1.0, 1.0));

    // Solid cast.
    assert_eq!(
        cuboid
            .cast_ray(&Isometry2::identity(), &ray_inside, std::f32::MAX, true)
            .unwrap(),
        0.0
    );

    // Non-solid cast.
    assert_eq!(
        cuboid
            .cast_ray(&Isometry2::identity(), &ray_inside, std::f32::MAX, false)
            .unwrap(),
        2.0
    );

    // The other ray does not intersect this shape.
    assert!(cuboid
        .cast_ray(&Isometry2::identity(), &ray_miss, std::f32::MAX, false)
        .is_none());
    assert!(cuboid
        .cast_ray(&Isometry2::identity(), &ray_miss, std::f32::MAX, true)
        .is_none());
}
