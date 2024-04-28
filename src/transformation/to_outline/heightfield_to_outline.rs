use crate::math::Real;
use crate::shape::{HeightField, HeightFieldCellStatus};
use na::Point3;

impl HeightField {
    /// Outlines this heightfield’s shape using polylines.
    pub fn to_outline(&self) -> (Vec<Point3<Real>>, Vec<[u32; 2]>) {
        todo!()
    }
}
