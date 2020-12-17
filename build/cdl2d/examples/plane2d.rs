extern crate nalgebra as na;

use cdl2d::shape::HalfSpace;
use na::Vector2;

fn main() {
    let halfspace = HalfSpace::new(Vector2::<f32>::y_axis());

    assert!(halfspace.normal.x == 0.0);
    assert!(halfspace.normal.y == 1.0);
}
