extern crate nalgebra as na;

use na::Vector3;
use parry3d::shape::HalfSpace;

fn main() {
    let halfspace = HalfSpace::new(Vector3::<f32>::y_axis());

    assert!(halfspace.normal.x == 0.0);
    assert!(halfspace.normal.y == 1.0);
    assert!(halfspace.normal.z == 0.0);
}
