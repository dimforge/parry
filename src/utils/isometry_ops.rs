use crate::math::{Isometry, Point, Real, SimdReal, Vector};
use crate::math::{IsometryT, VectorT};
use na::SimdComplexField;
use na::Unit; // for .abs()

#[cfg(not(feature = "std"))]
use na::ComplexField;

/// Extra operations with isometries.
pub trait IsometryOps<T> {
    /// Transform a vector by the absolute value of the homogeneous matrix
    /// equivalent to `self`.
    fn absolute_transform_vector(&self, v: &VectorT<T>) -> VectorT<T>;
}

impl IsometryOps<Real> for Isometry {
    #[inline]
    fn absolute_transform_vector(&self, v: &Vector) -> Vector {
        self.rotation.to_rotation_matrix().into_inner().abs() * *v
    }
}

impl IsometryOps<SimdReal> for IsometryT<SimdReal> {
    #[inline]
    fn absolute_transform_vector(&self, v: &VectorT<SimdReal>) -> VectorT<SimdReal> {
        self.rotation
            .to_rotation_matrix()
            .into_inner()
            .map(|e| e.simd_abs())
            * *v
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
    fn transform_unit_vector(self, v: &Unit<Vector>) -> Unit<Vector>;
    /// Computes `self.inverse() * p`.
    fn inverse_transform_point(self, p: &Point) -> Point;
    /// Computes `self.inverse() * v`.
    fn inverse_transform_vector(self, v: &Vector) -> Vector;
    /// Computes `self.inverse() * v`.
    fn inverse_transform_unit_vector(self, v: &Unit<Vector>) -> Unit<Vector>;
}

impl IsometryOpt for Option<&Isometry> {
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
            rhs * iso
        } else {
            *rhs
        }
    }

    #[inline]
    fn transform_point(self, p: &Point) -> Point {
        if let Some(iso) = self {
            iso * p
        } else {
            *p
        }
    }

    #[inline]
    fn transform_vector(self, v: &Vector) -> Vector {
        if let Some(iso) = self {
            iso * v
        } else {
            *v
        }
    }

    #[inline]
    fn transform_unit_vector(self, v: &Unit<Vector>) -> Unit<Vector> {
        if let Some(iso) = self {
            iso * v
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
    fn inverse_transform_unit_vector(self, v: &Unit<Vector>) -> Unit<Vector> {
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
            rhs * iso
        } else {
            *rhs
        }
    }

    #[inline]
    fn transform_point(self, p: &Point) -> Point {
        if let Some(iso) = self {
            iso * p
        } else {
            *p
        }
    }

    #[inline]
    fn transform_vector(self, v: &Vector) -> Vector {
        if let Some(iso) = self {
            iso * v
        } else {
            *v
        }
    }

    #[inline]
    fn transform_unit_vector(self, v: &Unit<Vector>) -> Unit<Vector> {
        if let Some(iso) = self {
            iso * v
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
    fn inverse_transform_unit_vector(self, v: &Unit<Vector>) -> Unit<Vector> {
        if let Some(iso) = self {
            iso.inverse_transform_unit_vector(v)
        } else {
            *v
        }
    }
}
