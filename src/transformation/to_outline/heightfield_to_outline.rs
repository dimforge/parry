use crate::math::Real;
use crate::shape::{GenericHeightField, HeightFieldCellStatus, HeightFieldStorage};
use na::Point3;

impl<Heights, Status> GenericHeightField<Heights, Status>
where
    Heights: HeightFieldStorage<Item = Real>,
    Status: HeightFieldStorage<Item = HeightFieldCellStatus>,
{
    /// Discretize the boundary of this ball as a triangle-mesh.
    pub fn to_outline(&self) -> (Vec<Point3<Real>>, Vec<[u32; 2]>) {
        todo!()
    }
}
