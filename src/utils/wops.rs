//! Miscellaneous utilities.

use crate::math::Real;
use crate::simd::{SimdBool, SimdReal};
use na::{Matrix3, Point2, Point3, Scalar, SimdRealField, Vector2, Vector3};
use simba::simd::SimdValue;

#[cfg(feature = "simd-is-enabled")]
use na::SimdPartialOrd;

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

pub(crate) trait WComponent: Sized {
    type Element;

    fn min_component(self) -> Self::Element;
    fn max_component(self) -> Self::Element;
}

impl WComponent for Real {
    type Element = Real;

    fn min_component(self) -> Self::Element {
        self
    }
    fn max_component(self) -> Self::Element {
        self
    }
}

#[cfg(feature = "simd-is-enabled")]
impl WComponent for SimdReal {
    type Element = Real;

    fn min_component(self) -> Self::Element {
        self.simd_horizontal_min()
    }
    fn max_component(self) -> Self::Element {
        self.simd_horizontal_max()
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

pub(crate) trait WVec: Sized {
    type Element;

    fn horizontal_inf(&self) -> Self::Element;
    fn horizontal_sup(&self) -> Self::Element;
}

impl<N: Scalar + Copy + WComponent> WVec for Vector2<N>
where
    N::Element: Scalar,
{
    type Element = Vector2<N::Element>;

    fn horizontal_inf(&self) -> Self::Element {
        Vector2::new(self.x.min_component(), self.y.min_component())
    }

    fn horizontal_sup(&self) -> Self::Element {
        Vector2::new(self.x.max_component(), self.y.max_component())
    }
}

impl<N: Scalar + Copy + WComponent> WVec for Point2<N>
where
    N::Element: Scalar,
{
    type Element = Point2<N::Element>;

    fn horizontal_inf(&self) -> Self::Element {
        Point2::new(self.x.min_component(), self.y.min_component())
    }

    fn horizontal_sup(&self) -> Self::Element {
        Point2::new(self.x.max_component(), self.y.max_component())
    }
}

impl<N: Scalar + Copy + WComponent> WVec for Vector3<N>
where
    N::Element: Scalar,
{
    type Element = Vector3<N::Element>;

    fn horizontal_inf(&self) -> Self::Element {
        Vector3::new(
            self.x.min_component(),
            self.y.min_component(),
            self.z.min_component(),
        )
    }

    fn horizontal_sup(&self) -> Self::Element {
        Vector3::new(
            self.x.max_component(),
            self.y.max_component(),
            self.z.max_component(),
        )
    }
}

impl<N: Scalar + Copy + WComponent> WVec for Point3<N>
where
    N::Element: Scalar,
{
    type Element = Point3<N::Element>;

    fn horizontal_inf(&self) -> Self::Element {
        Point3::new(
            self.x.min_component(),
            self.y.min_component(),
            self.z.min_component(),
        )
    }

    fn horizontal_sup(&self) -> Self::Element {
        Point3::new(
            self.x.max_component(),
            self.y.max_component(),
            self.z.max_component(),
        )
    }
}

pub(crate) trait WCrossMatrix: Sized {
    type CrossMat;

    fn gcross_matrix(self) -> Self::CrossMat;
}

impl WCrossMatrix for Vector3<Real> {
    type CrossMat = Matrix3<Real>;

    #[inline]
    #[rustfmt::skip]
    fn gcross_matrix(self) -> Self::CrossMat {
        Matrix3::new(
            0.0, -self.z, self.y,
            self.z, 0.0, -self.x,
            -self.y, self.x, 0.0,
        )
    }
}

impl WCrossMatrix for Vector2<Real> {
    type CrossMat = Vector2<Real>;

    #[inline]
    fn gcross_matrix(self) -> Self::CrossMat {
        Vector2::new(-self.y, self.x)
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

pub(crate) trait WDot<Rhs>: Sized {
    type Result;
    fn gdot(&self, rhs: Rhs) -> Self::Result;
}

impl WDot<Vector3<Real>> for Vector3<Real> {
    type Result = Real;

    fn gdot(&self, rhs: Vector3<Real>) -> Self::Result {
        self.x * rhs.x + self.y * rhs.y + self.z * rhs.z
    }
}

impl WDot<Vector2<Real>> for Vector2<Real> {
    type Result = Real;

    fn gdot(&self, rhs: Vector2<Real>) -> Self::Result {
        self.x * rhs.x + self.y * rhs.y
    }
}

impl WDot<Real> for Real {
    type Result = Real;

    fn gdot(&self, rhs: Real) -> Self::Result {
        *self * rhs
    }
}

#[cfg(feature = "simd-is-enabled")]
impl WCrossMatrix for Vector3<SimdReal> {
    type CrossMat = Matrix3<SimdReal>;

    #[inline]
    #[rustfmt::skip]
    fn gcross_matrix(self) -> Self::CrossMat {
        Matrix3::new(
            num::zero(), -self.z, self.y,
            self.z, num::zero(), -self.x,
            -self.y, self.x, num::zero(),
        )
    }
}

#[cfg(feature = "simd-is-enabled")]
impl WCrossMatrix for Vector2<SimdReal> {
    type CrossMat = Vector2<SimdReal>;

    #[inline]
    fn gcross_matrix(self) -> Self::CrossMat {
        Vector2::new(-self.y, self.x)
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

#[cfg(feature = "simd-is-enabled")]
impl WDot<Vector3<SimdReal>> for Vector3<SimdReal> {
    type Result = SimdReal;

    fn gdot(&self, rhs: Vector3<SimdReal>) -> Self::Result {
        self.x * rhs.x + self.y * rhs.y + self.z * rhs.z
    }
}

#[cfg(feature = "simd-is-enabled")]
impl WDot<Vector2<SimdReal>> for Vector2<SimdReal> {
    type Result = SimdReal;

    fn gdot(&self, rhs: Vector2<SimdReal>) -> Self::Result {
        self.x * rhs.x + self.y * rhs.y
    }
}

#[cfg(feature = "simd-is-enabled")]
impl WDot<SimdReal> for SimdReal {
    type Result = SimdReal;

    fn gdot(&self, rhs: SimdReal) -> Self::Result {
        *self * rhs
    }
}
