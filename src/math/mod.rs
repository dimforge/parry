//! Linear algebra type aliases.

#[cfg(feature = "linalg-glam")]
mod glam_impls_2d;
#[cfg(feature = "linalg-glam")]
mod glam_impls_simd_2d;
#[cfg(all(feature = "dim2", feature = "linalg-glam"))]
mod math_glam_2d;

#[cfg(all(feature = "dim2", feature = "linalg-glam"))]
pub use math_glam_2d::*;
#[cfg(feature = "linalg-glam")]
pub use {glam_impls_2d::*, glam_impls_simd_2d::*};

#[cfg(feature = "linalg-glam")]
mod glam_impls_3d;
#[cfg(feature = "linalg-glam")]
mod glam_impls_simd_3d;
#[cfg(all(feature = "dim3", feature = "linalg-glam"))]
mod math_glam_3d;

#[cfg(feature = "linalg-glam")]
pub use {glam_impls_3d::*, glam_impls_simd_3d::*};

#[cfg(all(feature = "dim3", feature = "linalg-glam"))]
pub use math_glam_3d::*;

use std::borrow::Borrow;

#[cfg(all(feature = "dim2", feature = "linalg-nalgebra"))]
mod math_nalgebra_2d;
#[cfg(all(feature = "dim2", feature = "linalg-nalgebra"))]
pub use math_nalgebra_2d::*;

#[cfg(all(feature = "dim3", feature = "linalg-nalgebra"))]
mod math_nalgebra_3d;

mod isometry_ops;

pub(crate) use crate::utils::{WBasis, WCross, WSign};

pub use isometry_ops::*;
#[cfg(all(feature = "dim3", feature = "linalg-nalgebra"))]
pub use math_nalgebra_3d::*;

pub use simd::*;

/// The scalar type used throughout this crate.
#[cfg(feature = "f64")]
pub type Real = f64;

/// The scalar type used throughout this crate.
#[cfg(feature = "f32")]
pub type Real = f32;

#[cfg(not(feature = "simd-is-enabled"))]
mod simd {
    use simba::simd::AutoBoolx4;
    /// The number of lanes of a SIMD number.
    pub const SIMD_WIDTH: usize = 4;
    /// SIMD_WIDTH - 1
    pub const SIMD_LAST_INDEX: usize = 3;

    /// A SIMD float with SIMD_WIDTH lanes.
    #[cfg(feature = "f32")]
    pub type SimdReal = simba::simd::AutoF32x4;

    /// A SIMD float with SIMD_WIDTH lanes.
    #[cfg(feature = "f64")]
    pub type SimdReal = simba::simd::AutoF64x4;

    /// A SIMD bool with SIMD_WIDTH lanes.
    pub type SimdBool = AutoBoolx4;
}

#[cfg(feature = "simd-is-enabled")]
mod simd {
    #[cfg(all(feature = "simd-nightly", feature = "f32"))]
    pub use simba::simd::{f32x4 as SimdReal, m32x4 as SimdBool};
    #[cfg(all(feature = "simd-stable", feature = "f32"))]
    pub use simba::simd::{WideBoolF32x4 as SimdBool, WideF32x4 as SimdReal};

    #[cfg(all(feature = "simd-nightly", feature = "f64"))]
    pub use simba::simd::{f64x4 as SimdReal, m64x4 as SimdBool};
    #[cfg(all(feature = "simd-stable", feature = "f64"))]
    pub use simba::simd::{WideBoolF64x4 as SimdBool, WideF64x4 as SimdReal};

    /// The number of lanes of a SIMD number.
    pub const SIMD_WIDTH: usize = 4;
    /// SIMD_WIDTH - 1
    pub const SIMD_LAST_INDEX: usize = 3;
}

#[cfg(feature = "linalg-glam")]
pub trait IntoInner: Sized {
    #[inline]
    fn into_inner(self) -> Self {
        self
    }
}

#[cfg(feature = "linalg-glam")]
pub trait Mat2Ops {
    type Element;

    /// Assuming Self represents a 2D rotation, compute the rotation angle.
    fn angle(&self) -> Self::Element;
}

#[cfg(feature = "linalg-glam")]
pub trait RotationOps: Sized {
    type Vector;
    type UnitVector;
    type AngVector;
    type RotationMatrix;
    type RotationBetween;

    fn identity() -> Self;
    fn from_scaled_axis(axis_angle: Self::AngVector) -> Self;
    fn new_eps(axis_angle: Self::AngVector, val: Real) -> Self;
    fn from_rotation_matrix(mat: &Self::RotationMatrix) -> Self;
    fn to_rotation_matrix(&self) -> Self::RotationMatrix;
    fn as_mut_unchecked(&mut self) -> &mut Self {
        self
    }
    fn renormalize(&mut self);
    fn renormalize_fast(&mut self);
    fn rotation_between(a: &Self::Vector, b: &Self::Vector) -> Self::RotationBetween;
    fn scaled_rotation_between_axis(
        a: &Self::Vector,
        b: &Self::Vector,
        s: f32,
    ) -> Self::RotationBetween;
    fn scaled_axis(&self) -> Self::AngVector;
    fn axis_angle(&self) -> Option<(Self::UnitVector, Real)>;
    fn imag(&self) -> Self::AngVector;
}

#[cfg(feature = "linalg-glam")]
pub trait PointOps: Sized {
    type Coords;
    // Not exactly the recommended signature for an `as_` function, but does the trick
    // to make it work with glam.
    #[inline(always)]
    fn as_vector(self) -> Self::Coords {
        self.into_vector()
    }
    fn into_vector(self) -> Self::Coords;
    fn as_vector_mut(&mut self) -> &mut Self::Coords;
    fn origin() -> Self;
}

#[cfg(feature = "linalg-glam")]
pub trait VectorU32Ops {
    type Vector;

    fn zeros() -> Self;
    fn repeat(val: u32) -> Self;
    fn inf(&self, rhs: &Self) -> Self;
    fn sup(&self, rhs: &Self) -> Self;
    fn map(&self, f: impl FnMut(u32) -> Real) -> Self::Vector;
    fn apply(&mut self, f: impl FnMut(&mut u32));
}

pub trait VectorOps: Sized {
    fn try_normalize_eps(self, eps: Real) -> Option<Self>;
}

impl VectorOps for na::Vector3<Real> {
    fn try_normalize_eps(self, eps: Real) -> Option<Self> {
        self.try_normalize(eps)
    }
}

impl VectorOps for na::Vector2<Real> {
    fn try_normalize_eps(self, eps: Real) -> Option<Self> {
        self.try_normalize(eps)
    }
}

#[cfg(feature = "linalg-glam")] // To define a common interface between glam and nalgebra.
pub trait MatrixOps<const DIM: usize>: Sized {
    type Vector;
    type SymmetricEigen;

    fn zeros() -> Self;
    fn from_columns(columns: &[Self::Vector; DIM]) -> Self;
    fn column_mut(&mut self, col: usize) -> &mut Self::Vector;
    fn from_matrix_unchecked(mat: Self) -> Self {
        mat
    }

    fn from_diagonal_element(elt: Real) -> Self;

    #[cfg(feature = "dim2")]
    fn new(m00: f32, m01: f32, m10: f32, m11: f32) -> Self;
    #[cfg(feature = "dim3")]
    fn new(
        m00: f32,
        m01: f32,
        m02: f32,
        m10: f32,
        m11: f32,
        m12: f32,
        m20: f32,
        m21: f32,
        m22: f32,
    ) -> Self;
    fn try_inverse(&self) -> Option<Self>;
    fn abs(&self) -> Self;
    fn symmetric_eigen(&self) -> Self::SymmetricEigen;
    fn swap_columns(&mut self, i: usize, j: usize);
    fn column(&self, i: usize) -> Self::Vector;
    fn symmetric_eigenvalues(&self) -> Self::Vector;
}

#[cfg(feature = "linalg-glam")] // To define a common interface between glam and nalgebra.
pub trait GlamVectorOps<const DIM: usize>: Sized {
    type Matrix;

    fn x() -> Self;
    fn y() -> Self;
    #[cfg(feature = "dim3")]
    fn z() -> Self;
    fn x_axis() -> Self {
        Self::x()
    }
    fn y_axis() -> Self {
        Self::y()
    }
    #[cfg(feature = "dim3")]
    fn z_axis() -> Self {
        Self::z()
    }
    fn ith(i: usize, val: Real) -> Self;
    fn ith_axis(i: usize) -> Self {
        Self::ith(i, 1.0)
    }
    fn from_row_slice(elts: &[Real]) -> Self;
    fn repeat(value: Real) -> Self;
    fn orthonormal_subspace_basis(dirs: &[Self], f: impl FnMut(&Self) -> bool);
    fn iter(&self) -> std::array::IntoIter<&Real, DIM>;
    #[inline(always)]
    fn cast<T>(self) -> Self {
        self
    }
    #[inline(always)]
    fn select(self, take_self: bool, other: Self) -> Self {
        if take_self {
            self
        } else {
            other
        }
    }
    fn as_slice(&self) -> &[Real];
    fn map(&self, f: impl FnMut(Real) -> Real) -> Self;
    fn zip_apply<F>(&mut self, b: &Self, f: F)
    where
        F: FnMut(&mut Real, Real);
    fn norm_squared(&self) -> Real;
    fn norm(&self) -> Real;
    fn inf(&self, rhs: &Self) -> Self;
    fn sup(&self, rhs: &Self) -> Self;
    fn component_div(&self, rhs: &Self) -> Self;
    fn component_mul(&self, rhs: &Self) -> Self;
    fn component_mul_assign(&mut self, rhs: &Self);
    fn zeros() -> Self;
    fn is_zero(self) -> bool;
    #[inline(always)]
    fn into_owned(self) -> Self {
        self
    }
    #[inline(always)]
    fn normalize(mut self) -> Self {
        let _ = self.normalize_mut();
        self
    }
    fn normalize_mut(&mut self) -> Real;
    fn try_normalize_mut(&mut self, eps: Real) -> Option<Real>;
    fn fill(&mut self, val: Real) {
        *self = Self::repeat(val)
    }
    fn imin(&self) -> usize;
    fn iamin(&self) -> usize;
    fn iamax(&self) -> usize;
    fn amin(&self) -> Real;
    fn amax(&self) -> Real;

    fn outer_product(&self, rhs: &Self) -> Self::Matrix;
    fn angle(&self, rhs: &Self) -> Real;
    fn swap_rows(&mut self, i: usize, j: usize);

    /*
     * UnitVector
     */
    fn try_new_and_get(v: Self, eps: Real) -> Option<(Self, Real)>;

    #[inline]
    fn try_new(v: Self, eps: Real) -> Option<Self> {
        Self::try_new_and_get(v, eps).map(|res| res.0)
    }

    #[inline]
    fn new_normalize(v: Self) -> Self {
        v.normalize()
    }

    #[inline]
    fn new_unchecked(v: Self) -> Self {
        v
    }
}

#[inline]
pub fn distance(a: impl Borrow<Point>, b: impl Borrow<Point>) -> Real {
    (*a.borrow() - *b.borrow()).norm()
}

#[inline]
pub fn distance_squared(a: impl Borrow<Point>, b: impl Borrow<Point>) -> Real {
    (*a.borrow() - *b.borrow()).norm_squared()
}

#[inline]
pub fn center(a: impl Borrow<Point>, b: impl Borrow<Point>) -> Point {
    (*a.borrow() + b.borrow().as_vector()) / 2.0
}

#[cfg(feature = "linalg-nalgebra")]
#[inline]
pub fn outer_product(a: impl Borrow<Vector>, b: impl Borrow<Vector>) -> Matrix {
    a.borrow() * b.borrow().transpose()
}

#[cfg(feature = "linalg-glam")]
#[inline]
pub fn outer_product(a: impl Borrow<Vector>, b: impl Borrow<Vector>) -> Matrix {
    a.borrow().outer_product(b.borrow())
}

#[cfg(feature = "linalg-glam")]
#[derive(Copy, Clone, PartialEq, Debug)]
pub struct SymmetricEigen2 {
    pub eigenvectors: glam::Mat2,
    pub eigenvalues: glam::Vec2,
}

#[cfg(feature = "linalg-glam")]
#[derive(Copy, Clone, PartialEq, Debug)]
pub struct SymmetricEigen3 {
    pub eigenvectors: glam::Mat3,
    pub eigenvalues: glam::Vec3,
}
