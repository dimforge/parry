use parry2d::math::Vector;
use parry2d::shape::Cuboid;

fn main() {
    let cuboid = Cuboid::new(Vector::new(2.0, 1.0));

    assert!(cuboid.half_extents.x == 2.0);
    assert!(cuboid.half_extents.y == 1.0);
}
