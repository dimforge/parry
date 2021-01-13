use crate::math::Real;
use crate::shape::Tetrahedron;
use bitflags::_core::cmp::min;
use na::{Point3, Vector3};
use std::convert::TryFrom;

#[derive(Copy, Clone, Debug)]
pub enum Axis {
    X = 0,
    Y,
    Z,
}

impl Axis {
    pub fn try_from(val: usize) -> Option<Self> {
        match val {
            0 => Some(Axis::X),
            1 => Some(Axis::Y),
            2 => Some(Axis::Z),
            _ => None,
        }
    }
}

#[derive(Copy, Clone, Debug)]
pub struct Plane {
    pub abc: Vector3<Real>,
    pub d: Real,
    pub axis: Axis,
    pub index: u32,
}
