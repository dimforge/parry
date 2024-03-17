use super::{
    Isometry, Point, Rotation, SimdIsometry, SimdRotation, SimdVector, Translation, UnitVector,
    Vector,
};
use na::SimdComplexField; // for .abs()

#[cfg(not(feature = "std"))]
use na::ComplexField;

#[cfg(feature = "linalg-glam")]
use super::MatrixOps;

/// Extra operations with isometries.
pub trait IsometryOps<V, R> {
    /// Transform a vector by the absolute value of the homogeneous matrix
    /// equivalent to `self`.
    fn absolute_transform_vector(&self, v: &V) -> V;
    fn rotation(&self) -> &R;
}

#[cfg(feature = "linalg-nalgebra")]
impl IsometryOps<Vector, Rotation> for Isometry {
    #[inline]
    fn absolute_transform_vector(&self, v: &Vector) -> Vector {
        self.rotation.to_rotation_matrix().into_inner().abs() * *v
    }

    #[inline(always)]
    fn rotation(&self) -> &Rotation {
        &self.rotation
    }
}

#[cfg(feature = "linalg-glam")]
impl IsometryOps<Vector, Rotation> for Isometry {
    #[cfg(feature = "dim2")]
    #[inline]
    fn absolute_transform_vector(&self, v: &Vector) -> Vector {
        self.rotation.abs() * *v
    }

    #[cfg(feature = "dim3")]
    #[inline]
    fn absolute_transform_vector(&self, v: &Vector) -> Vector {
        super::RotationMatrix::from_quat(self.rotation).abs() * *v
    }

    #[inline(always)]
    fn rotation(&self) -> &Rotation {
        &self.rotation
    }
}

impl IsometryOps<SimdVector, SimdRotation> for SimdIsometry {
    #[inline]
    fn absolute_transform_vector(&self, v: &SimdVector) -> SimdVector {
        self.rotation
            .to_rotation_matrix()
            .into_inner()
            .map(|e| e.simd_abs())
            * *v
    }

    #[inline(always)]
    fn rotation(&self) -> &SimdRotation {
        &self.rotation
    }
}

/// Various operations usable with `Option<Isometry>` and `Option<&Isometry>`
/// where `None` is assumed to be equivalent to the identity.
pub trait IsometryOpt {
    /// Computes `self.inverse() * rhs`.
    fn inv_mul(self, rhs: &Isometry) -> Isometry;
    /// Computes `rhs * self`.
    fn prepend_to(self, rhs: &Isometry) -> Isometry;
    /// Computes `self * p`.
    fn transform_point(self, p: &Point) -> Point;
    /// Computes `self * v`.
    fn transform_vector(self, v: &Vector) -> Vector;
    /// Computes `self * v`.
    fn transform_unit_vector(self, v: &UnitVector) -> UnitVector;
    /// Computes `self.inverse() * p`.
    fn inverse_transform_point(self, p: &Point) -> Point;
    /// Computes `self.inverse() * v`.
    fn inverse_transform_vector(self, v: &Vector) -> Vector;
    /// Computes `self.inverse() * v`.
    fn inverse_transform_unit_vector(self, v: &UnitVector) -> UnitVector;
}

#[cfg(feature = "linalg-glam")]
impl<'a> IsometryOpt for Isometry {
    fn inv_mul(self, rhs: &Isometry) -> Isometry {
        self.inverse() * *rhs
    }
    fn prepend_to(self, rhs: &Isometry) -> Isometry {
        *rhs * self
    }
    fn transform_point(self, p: &Point) -> Point {
        self.rotation * *p + self.translation
    }
    fn transform_vector(self, v: &Vector) -> Vector {
        self.rotation * *v
    }
    fn transform_unit_vector(self, v: &UnitVector) -> UnitVector {
        self.rotation * *v
    }
    fn inverse_transform_point(self, p: &Point) -> Point {
        self.inverse().transform_point(p)
    }
    fn inverse_transform_vector(self, v: &Vector) -> Vector {
        self.rotation.inverse() * *v
    }
    fn inverse_transform_unit_vector(self, v: &UnitVector) -> UnitVector {
        self.rotation.inverse() * *v
    }
}

impl<'a> IsometryOpt for Option<&'a Isometry> {
    #[inline]
    fn inv_mul(self, rhs: &Isometry) -> Isometry {
        if let Some(iso) = self {
            iso.inv_mul(rhs)
        } else {
            *rhs
        }
    }

    #[inline]
    fn prepend_to(self, rhs: &Isometry) -> Isometry {
        if let Some(iso) = self {
            *rhs * *iso
        } else {
            *rhs
        }
    }

    #[inline]
    fn transform_point(self, p: &Point) -> Point {
        if let Some(iso) = self {
            iso.transform_point(p)
        } else {
            *p
        }
    }

    #[inline]
    fn transform_vector(self, v: &Vector) -> Vector {
        if let Some(iso) = self {
            iso.rotation * *v
        } else {
            *v
        }
    }

    #[inline]
    fn transform_unit_vector(self, v: &UnitVector) -> UnitVector {
        if let Some(iso) = self {
            iso.rotation * *v
        } else {
            *v
        }
    }

    #[inline]
    fn inverse_transform_point(self, p: &Point) -> Point {
        if let Some(iso) = self {
            iso.inverse_transform_point(p)
        } else {
            *p
        }
    }

    #[inline]
    fn inverse_transform_vector(self, v: &Vector) -> Vector {
        if let Some(iso) = self {
            iso.inverse_transform_vector(v)
        } else {
            *v
        }
    }

    #[inline]
    fn inverse_transform_unit_vector(self, v: &UnitVector) -> UnitVector {
        if let Some(iso) = self {
            iso.inverse_transform_unit_vector(v)
        } else {
            *v
        }
    }
}

impl IsometryOpt for Option<Isometry> {
    #[inline]
    fn inv_mul(self, rhs: &Isometry) -> Isometry {
        if let Some(iso) = self {
            iso.inv_mul(rhs)
        } else {
            *rhs
        }
    }

    #[inline]
    fn prepend_to(self, rhs: &Isometry) -> Isometry {
        if let Some(iso) = self {
            *rhs * iso
        } else {
            *rhs
        }
    }

    #[inline]
    fn transform_point(self, p: &Point) -> Point {
        if let Some(iso) = self {
            iso.transform_point(p)
        } else {
            *p
        }
    }

    #[inline]
    fn transform_vector(self, v: &Vector) -> Vector {
        if let Some(iso) = self {
            iso.rotation * *v
        } else {
            *v
        }
    }

    #[inline]
    fn transform_unit_vector(self, v: &UnitVector) -> UnitVector {
        if let Some(iso) = self {
            iso.rotation * *v
        } else {
            *v
        }
    }

    #[inline]
    fn inverse_transform_point(self, p: &Point) -> Point {
        if let Some(iso) = self {
            iso.inverse_transform_point(p)
        } else {
            *p
        }
    }

    #[inline]
    fn inverse_transform_vector(self, v: &Vector) -> Vector {
        if let Some(iso) = self {
            iso.inverse_transform_vector(v)
        } else {
            *v
        }
    }

    #[inline]
    fn inverse_transform_unit_vector(self, v: &UnitVector) -> UnitVector {
        if let Some(iso) = self {
            iso.inverse_transform_unit_vector(v)
        } else {
            *v
        }
    }
}
