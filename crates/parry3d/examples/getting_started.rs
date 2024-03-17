extern crate nalgebra as na;

use parry3d::math::{GlamVectorOps, Isometry, Point, Vector};
use parry3d::query::{Ray, RayCast};
use parry3d::shape::Cuboid;

fn main() {
    let cube = Cuboid::new(Vector::new(1.0f32, 1.0, 1.0));
    let ray = Ray::new(Point::new(0.0f32, 0.0, -1.0), Vector::z());

    assert!(cube.intersects_ray(&Isometry::identity(), &ray, std::f32::MAX));
}
