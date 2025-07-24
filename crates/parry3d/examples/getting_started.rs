extern crate nalgebra as na;

use na::{Isometry3, Point3, Vector3};
use parry3d::query::{Ray, RayCast};
use parry3d::shape::Cuboid;

fn main() {
    let cube = Cuboid::new(Vector3::new(1.0f32, 1.0, 1.0));
    let ray = Ray::new(Point3::new(0.0f32, 0.0, -1.0), Vector3::z());

    assert!(cube.intersects_ray(&Isometry3::identity(), &ray, f32::MAX, &()));
}
