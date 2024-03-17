use parry3d::math::{GlamVectorOps, Vector};
use parry3d::shape::HalfSpace;

fn main() {
    let halfspace = HalfSpace::new(Vector::y_axis());

    assert!(halfspace.normal.x == 0.0);
    assert!(halfspace.normal.y == 1.0);
    assert!(halfspace.normal.z == 0.0);
}
