use crate::math::{Point, Real, Vector};
use na::{self, Unit};
use std::mem;

/// Geometric description of a contact.
#[derive(Debug, PartialEq, Copy, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Contact {
    /// Position of the contact on the first object. The position is expressed in local-space of the first shape.
    pub point1: Point<Real>,

    /// Position of the contact on the second object. The position is expressed in local-space of the second shape.
    pub point2: Point<Real>,

    /// Contact normal expressed in the local-space of the first shape.
    ///
    /// This is an outward normal, i.e., it points towards the exterior of the first shape.
    pub normal1: Unit<Vector<Real>>,

    /// Contact normal expressed in the local-space of the second shape.
    ///
    /// This is an outward normal, i.e., it points towards the exterior of the second shape.
    pub normal2: Unit<Vector<Real>>,

    /// Distance between the two contact points.
    ///
    /// If this is negative, this contact represents a penetration.
    pub dist: Real,
}

impl Contact {
    /// Creates a new contact.
    #[inline]
    pub fn new(
        point1: Point<Real>,
        point2: Point<Real>,
        normal1: Unit<Vector<Real>>,
        normal2: Unit<Vector<Real>>,
        dist: Real,
    ) -> Self {
        Contact {
            point1,
            point2,
            normal1,
            normal2,
            dist,
        }
    }
}

impl Contact {
    /// Swaps the points and normals of this contact.
    #[inline]
    pub fn flip(&mut self) {
        mem::swap(&mut self.point1, &mut self.point2);
        mem::swap(&mut self.point1, &mut self.point2);
    }

    /// Returns a new contact containing the swapped points and normals of `self`.
    #[inline]
    pub fn flipped(mut self) -> Self {
        self.flip();
        self
    }
}
