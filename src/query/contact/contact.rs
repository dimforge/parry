use crate::math::{Isometry, Point, Real, Vector};
use na::{self, Unit};
use std::mem;

#[cfg(feature = "rkyv")]
use rkyv::{bytecheck, CheckBytes};

/// Geometric description of a contact.
#[derive(Debug, PartialEq, Copy, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Deserialize, rkyv::Serialize, CheckBytes),
    archive(as = "Self")
)]
pub struct Contact {
    /// Position of the contact on the first object.
    pub point1: Point,

    /// Position of the contact on the second object.
    pub point2: Point,

    /// Contact normal, pointing towards the exterior of the first shape.
    pub normal1: Unit<Vector>,

    /// Contact normal, pointing towards the exterior of the second shape.
    ///
    /// If these contact data are expressed in world-space, this normal is equal to `-normal1`.
    pub normal2: Unit<Vector>,

    /// Distance between the two contact points.
    ///
    /// If this is negative, this contact represents a penetration.
    pub dist: Real,
}

impl Contact {
    /// Creates a new contact.
    #[inline]
    pub fn new(
        point1: Point,
        point2: Point,
        normal1: Unit<Vector>,
        normal2: Unit<Vector>,
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
        mem::swap(&mut self.normal1, &mut self.normal2);
    }

    /// Returns a new contact containing the swapped points and normals of `self`.
    #[inline]
    pub fn flipped(mut self) -> Self {
        self.flip();
        self
    }

    /// Transform the points and normals from this contact by
    /// the given transformations.
    #[inline]
    pub fn transform_by_mut(&mut self, pos1: &Isometry, pos2: &Isometry) {
        self.point1 = pos1 * self.point1;
        self.point2 = pos2 * self.point2;
        self.normal1 = pos1 * self.normal1;
        self.normal2 = pos2 * self.normal2;
    }

    /// Transform `self.point1` and `self.normal1` by the `pos`.
    pub fn transform1_by_mut(&mut self, pos: &Isometry) {
        self.point1 = pos * self.point1;
        self.normal1 = pos * self.normal1;
    }
}
