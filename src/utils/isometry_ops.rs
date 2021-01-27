use crate::math::{Isometry, Point, Real, Vector};
use na::Unit;

/// Extra operations with isometries.
pub trait IsometryOps {
    /// Transform a vector by the absolute value of the homogeneous matrix
    /// equivalent to `self`.
    fn absolute_transform_vector(&self, v: &Vector<Real>) -> Vector<Real>;
}

impl IsometryOps for Isometry<Real> {
    #[inline]
    fn absolute_transform_vector(&self, v: &Vector<Real>) -> Vector<Real> {
        self.rotation.to_rotation_matrix().into_inner().abs() * *v
    }
}

/// Various operations usable with `Option<Isometry>` and `Option<&Isometry>`
/// where `None` is assumed to be equivalent to the identity.
pub trait IsometryOpt {
    /// Computes `self.inverse() * rhs`.
    fn inv_mul(self, rhs: &Isometry<Real>) -> Isometry<Real>;
    /// Computes `rhs * self`.
    fn prepend_to(self, rhs: &Isometry<Real>) -> Isometry<Real>;
    /// Computes `self * p`.
    fn transform_point(self, p: &Point<Real>) -> Point<Real>;
    /// Computes `self * v`.
    fn transform_vector(self, v: &Vector<Real>) -> Vector<Real>;
    /// Computes `self * v`.
    fn transform_unit_vector(self, v: &Unit<Vector<Real>>) -> Unit<Vector<Real>>;
    /// Computes `self.inverse() * p`.
    fn inverse_transform_point(self, p: &Point<Real>) -> Point<Real>;
    /// Computes `self.inverse() * v`.
    fn inverse_transform_vector(self, v: &Vector<Real>) -> Vector<Real>;
    /// Computes `self.inverse() * v`.
    fn inverse_transform_unit_vector(self, v: &Unit<Vector<Real>>) -> Unit<Vector<Real>>;
}

impl<'a> IsometryOpt for Option<&'a Isometry<Real>> {
    #[inline]
    fn inv_mul(self, rhs: &Isometry<Real>) -> Isometry<Real> {
        if let Some(iso) = self {
            iso.inv_mul(&rhs)
        } else {
            *rhs
        }
    }

    #[inline]
    fn prepend_to(self, rhs: &Isometry<Real>) -> Isometry<Real> {
        if let Some(iso) = self {
            rhs * iso
        } else {
            *rhs
        }
    }

    #[inline]
    fn transform_point(self, p: &Point<Real>) -> Point<Real> {
        if let Some(iso) = self {
            iso * p
        } else {
            *p
        }
    }

    #[inline]
    fn transform_vector(self, v: &Vector<Real>) -> Vector<Real> {
        if let Some(iso) = self {
            iso * v
        } else {
            *v
        }
    }

    #[inline]
    fn transform_unit_vector(self, v: &Unit<Vector<Real>>) -> Unit<Vector<Real>> {
        if let Some(iso) = self {
            iso * v
        } else {
            *v
        }
    }

    #[inline]
    fn inverse_transform_point(self, p: &Point<Real>) -> Point<Real> {
        if let Some(iso) = self {
            iso.inverse_transform_point(p)
        } else {
            *p
        }
    }

    #[inline]
    fn inverse_transform_vector(self, v: &Vector<Real>) -> Vector<Real> {
        if let Some(iso) = self {
            iso.inverse_transform_vector(v)
        } else {
            *v
        }
    }

    #[inline]
    fn inverse_transform_unit_vector(self, v: &Unit<Vector<Real>>) -> Unit<Vector<Real>> {
        if let Some(iso) = self {
            iso.inverse_transform_unit_vector(v)
        } else {
            *v
        }
    }
}

impl IsometryOpt for Option<Isometry<Real>> {
    #[inline]
    fn inv_mul(self, rhs: &Isometry<Real>) -> Isometry<Real> {
        if let Some(iso) = self {
            iso.inv_mul(&rhs)
        } else {
            *rhs
        }
    }

    #[inline]
    fn prepend_to(self, rhs: &Isometry<Real>) -> Isometry<Real> {
        if let Some(iso) = self {
            rhs * iso
        } else {
            *rhs
        }
    }

    #[inline]
    fn transform_point(self, p: &Point<Real>) -> Point<Real> {
        if let Some(iso) = self {
            iso * p
        } else {
            *p
        }
    }

    #[inline]
    fn transform_vector(self, v: &Vector<Real>) -> Vector<Real> {
        if let Some(iso) = self {
            iso * v
        } else {
            *v
        }
    }

    #[inline]
    fn transform_unit_vector(self, v: &Unit<Vector<Real>>) -> Unit<Vector<Real>> {
        if let Some(iso) = self {
            iso * v
        } else {
            *v
        }
    }

    #[inline]
    fn inverse_transform_point(self, p: &Point<Real>) -> Point<Real> {
        if let Some(iso) = self {
            iso.inverse_transform_point(p)
        } else {
            *p
        }
    }

    #[inline]
    fn inverse_transform_vector(self, v: &Vector<Real>) -> Vector<Real> {
        if let Some(iso) = self {
            iso.inverse_transform_vector(v)
        } else {
            *v
        }
    }

    #[inline]
    fn inverse_transform_unit_vector(self, v: &Unit<Vector<Real>>) -> Unit<Vector<Real>> {
        if let Some(iso) = self {
            iso.inverse_transform_unit_vector(v)
        } else {
            *v
        }
    }
}
