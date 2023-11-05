use crate::math::Real;
use crate::shape::{GenericHeightField, HeightFieldCellStatus, HeightFieldStorage};
use na::Point3;

#[cfg(feature = "alloc")]
use alloc::{vec::Vec};

impl<Storage: HeightFieldStorage> GenericHeightField<Storage> {
    /// Outlines this heightfield’s shape using polylines.
    pub fn to_outline(&self) -> (Vec<Point3<Real>>, Vec<[u32; 2]>) {
        todo!()
    }
}
