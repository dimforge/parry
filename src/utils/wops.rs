//! Miscellaneous utilities.

use crate::math::Real;
use crate::simd::{SimdBool, SimdReal};
use na::{Scalar, SimdRealField, Vector2, Vector3};
use simba::simd::SimdValue;

/// Conditionally swaps each lanes of `a` with those of `b`.
///
/// For each `i in [0..SIMD_WIDTH[`, if `do_swap.extract(i)` is `true` then
/// `a.extract(i)` is swapped with `b.extract(i)`.
pub fn simd_swap(do_swap: SimdBool, a: &mut SimdReal, b: &mut SimdReal) {
    let _a = *a;
    *a = b.select(do_swap, *a);
    *b = _a.select(do_swap, *b);
}

/// Trait to copy the sign of each component of one scalar/vector/matrix to another.
pub trait WSign<Rhs>: Sized {
    // See SIMD implementations of copy_sign there: https://stackoverflow.com/a/57872652
    /// Copy the sign of each component of `self` to the corresponding component of `to`.
    fn copy_sign_to(self, to: Rhs) -> Rhs;
}

impl WSign<Real> for Real {
    fn copy_sign_to(self, to: Self) -> Self {
        let minus_zero: Real = -0.0;
        let signbit = minus_zero.to_bits();
        Real::from_bits((signbit & self.to_bits()) | ((!signbit) & to.to_bits()))
    }
}

impl<N: Scalar + Copy + WSign<N>> WSign<Vector2<N>> for N {
    fn copy_sign_to(self, to: Vector2<N>) -> Vector2<N> {
        Vector2::new(self.copy_sign_to(to.x), self.copy_sign_to(to.y))
    }
}

impl<N: Scalar + Copy + WSign<N>> WSign<Vector3<N>> for N {
    fn copy_sign_to(self, to: Vector3<N>) -> Vector3<N> {
        Vector3::new(
            self.copy_sign_to(to.x),
            self.copy_sign_to(to.y),
            self.copy_sign_to(to.z),
        )
    }
}

impl<N: Scalar + Copy + WSign<N>> WSign<Vector2<N>> for Vector2<N> {
    fn copy_sign_to(self, to: Vector2<N>) -> Vector2<N> {
        Vector2::new(self.x.copy_sign_to(to.x), self.y.copy_sign_to(to.y))
    }
}

impl<N: Scalar + Copy + WSign<N>> WSign<Vector3<N>> for Vector3<N> {
    fn copy_sign_to(self, to: Vector3<N>) -> Vector3<N> {
        Vector3::new(
            self.x.copy_sign_to(to.x),
            self.y.copy_sign_to(to.y),
            self.z.copy_sign_to(to.z),
        )
    }
}

#[cfg(feature = "simd-is-enabled")]
impl WSign<SimdReal> for SimdReal {
    fn copy_sign_to(self, to: SimdReal) -> SimdReal {
        to.simd_copysign(self)
    }
}

/// Trait to compute the orthonormal basis of a vector.
pub trait WBasis: Sized {
    /// The type of the array of orthonormal vectors.
    type Basis;
    /// Computes the vectors which, when combined with `self`, form an orthonormal basis.
    fn orthonormal_basis(self) -> Self::Basis;
}

impl<N: SimdRealField + Copy> WBasis for Vector2<N> {
    type Basis = [Vector2<N>; 1];
    fn orthonormal_basis(self) -> [Vector2<N>; 1] {
        [Vector2::new(-self.y, self.x)]
    }
}

impl<N: SimdRealField + Copy + WSign<N>> WBasis for Vector3<N> {
    type Basis = [Vector3<N>; 2];
    // Robust and branchless implementation from Pixar:
    // https://graphics.pixar.com/library/OrthonormalB/paper.pdf
    fn orthonormal_basis(self) -> [Vector3<N>; 2] {
        let sign = self.z.copy_sign_to(N::one());
        let a = -N::one() / (sign + self.z);
        let b = self.x * self.y * a;

        [
            Vector3::new(
                N::one() + sign * self.x * self.x * a,
                sign * b,
                -sign * self.x,
            ),
            Vector3::new(b, sign + self.y * self.y * a, -self.y),
        ]
    }
}

pub(crate) trait WCross<Rhs>: Sized {
    type Result;
    fn gcross(&self, rhs: Rhs) -> Self::Result;
}

impl WCross<Vector3<Real>> for Vector3<Real> {
    type Result = Self;

    fn gcross(&self, rhs: Vector3<Real>) -> Self::Result {
        self.cross(&rhs)
    }
}

impl WCross<Vector2<Real>> for Vector2<Real> {
    type Result = Real;

    fn gcross(&self, rhs: Vector2<Real>) -> Self::Result {
        self.x * rhs.y - self.y * rhs.x
    }
}

impl WCross<Vector2<Real>> for Real {
    type Result = Vector2<Real>;

    fn gcross(&self, rhs: Vector2<Real>) -> Self::Result {
        Vector2::new(-rhs.y * *self, rhs.x * *self)
    }
}

#[cfg(feature = "simd-is-enabled")]
impl WCross<Vector3<SimdReal>> for Vector3<SimdReal> {
    type Result = Vector3<SimdReal>;

    fn gcross(&self, rhs: Self) -> Self::Result {
        self.cross(&rhs)
    }
}

#[cfg(feature = "simd-is-enabled")]
impl WCross<Vector2<SimdReal>> for SimdReal {
    type Result = Vector2<SimdReal>;

    fn gcross(&self, rhs: Vector2<SimdReal>) -> Self::Result {
        Vector2::new(-rhs.y * *self, rhs.x * *self)
    }
}

#[cfg(feature = "simd-is-enabled")]
impl WCross<Vector2<SimdReal>> for Vector2<SimdReal> {
    type Result = SimdReal;

    fn gcross(&self, rhs: Self) -> Self::Result {
        let yx = Vector2::new(rhs.y, rhs.x);
        let prod = self.component_mul(&yx);
        prod.x - prod.y
    }
}
