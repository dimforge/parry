extern crate nalgebra as na;

use na::Vector2;
use parry2d::shape::HalfSpace;

fn main() {
    let halfspace = HalfSpace::new(Vector2::<f32>::y_axis());

    assert!(halfspace.normal.x == 0.0);
    assert!(halfspace.normal.y == 1.0);
}
