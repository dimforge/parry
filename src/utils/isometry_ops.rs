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

pub trait IsometryOpt {
    fn inv_mul(self, rhs: &Isometry<Real>) -> Isometry<Real>;
    fn transform_point(self, p: &Point<Real>) -> Point<Real>;
    fn transform_vector(self, p: &Vector<Real>) -> Vector<Real>;
    fn transform_unit_vector(self, p: &Unit<Vector<Real>>) -> Unit<Vector<Real>>;
    fn inverse_transform_point(self, p: &Point<Real>) -> Point<Real>;
    fn inverse_transform_vector(self, v: &Vector<Real>) -> Vector<Real>;
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
