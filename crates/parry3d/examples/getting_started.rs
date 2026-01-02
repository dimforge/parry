use parry3d::math::{Pose, Vector};
use parry3d::query::{Ray, RayCast};
use parry3d::shape::Cuboid;

fn main() {
    let cube = Cuboid::new(Vector::new(1.0f32, 1.0, 1.0));
    let ray = Ray::new(Vector::new(0.0f32, 0.0, -1.0), Vector::Z);

    assert!(cube.intersects_ray(&Pose::identity(), &ray, f32::MAX));
}
