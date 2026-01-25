use parry2d::math::Vector;
use parry2d::shape::HalfSpace;

fn main() {
    let halfspace = HalfSpace::new(Vector::Y);

    assert!(halfspace.normal.x == 0.0);
    assert!(halfspace.normal.y == 1.0);
}
