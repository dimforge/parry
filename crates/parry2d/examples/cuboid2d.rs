extern crate nalgebra as na;

use na::Vector2;
use parry2d::shape::Cuboid;

fn main() {
    let cuboid = Cuboid::new(Vector2::new(2.0, 1.0));

    assert!(cuboid.half_extents.x == 2.0);
    assert!(cuboid.half_extents.y == 1.0);
}
