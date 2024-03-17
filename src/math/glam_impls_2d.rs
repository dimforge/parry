pub use super::{
    GlamVectorOps, MatrixOps, PointOps, Real, RotationMatrix, RotationOps, SimdIso2, SimdMat2,
    SymmetricEigen2, VectorOps, VectorU32Ops,
};
use crate::math::{IntoInner, Mat2Ops, WCross};
use crate::utils::{SdpMatrix2, WBasis, WSign};
use core::ops::Mul;
use glam::{Mat2, UVec2, Vec2};
use na::Vector2;
use simba::simd::SimdRealField;

/// A 2D isometry.
#[derive(Copy, Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde-serialize", derive(Serialize, Deserialize))]
pub struct Iso2 {
    pub rotation: Mat2,
    pub translation: Vec2,
}

impl Default for Iso2 {
    #[inline]
    fn default() -> Self {
        Self::identity()
    }
}

impl Iso2 {
    #[inline]
    pub fn new(translation: Vec2, axis_angle: Real) -> Self {
        Self {
            rotation: Mat2::from_angle(axis_angle),
            translation,
        }
    }

    #[inline]
    pub fn from_parts(translation: Vec2, rotation: Mat2) -> Self {
        Self {
            rotation,
            translation,
        }
    }

    #[inline]
    pub fn identity() -> Self {
        Self {
            rotation: Mat2::IDENTITY,
            translation: Vec2::ZERO,
        }
    }

    #[inline]
    pub fn rotation(axis_angle: Real) -> Self {
        Self {
            rotation: Mat2::from_angle(axis_angle),
            translation: Vec2::ZERO,
        }
    }

    #[inline]
    pub fn translation(x: Real, y: Real) -> Self {
        Self {
            rotation: Mat2::IDENTITY,
            translation: Vec2::new(x, y),
        }
    }

    #[inline]
    pub fn inverse(&self) -> Self {
        let inv_rot = self.rotation.inverse();
        Self {
            rotation: inv_rot,
            translation: inv_rot * (-self.translation),
        }
    }

    #[inline]
    pub fn inv_mul(&self, rhs: &Self) -> Self {
        self.inverse() * *rhs
    }

    #[inline]
    pub fn lerp_slerp(&self, rhs: Self, t: Real) -> Self {
        let delta = rhs.rotation * self.rotation.transpose();

        Self {
            rotation: self.rotation * Mat2::from_angle(delta.angle() * t),
            translation: self.translation.lerp(rhs.translation, t),
        }
    }
}

impl Mat2Ops for Mat2 {
    type Element = Real;
    #[inline]
    fn angle(&self) -> Real {
        use na::SimdRealField;
        let arr = self.as_ref();
        arr[1].clone().simd_atan2(arr[0].clone())
    }
}

impl From<Vec2> for Iso2 {
    #[inline]
    fn from(translation: Vec2) -> Self {
        Self {
            rotation: Mat2::IDENTITY,
            translation,
        }
    }
}

impl Mul<Iso2> for Mat2 {
    type Output = Iso2;
    #[inline]
    fn mul(self, rhs: Iso2) -> Iso2 {
        Iso2 {
            rotation: self * rhs.rotation,
            translation: self * rhs.translation,
        }
    }
}

impl Mul<Iso2> for Iso2 {
    type Output = Iso2;
    #[inline]
    fn mul(self, rhs: Iso2) -> Iso2 {
        &self * &rhs
    }
}

impl<'a> Mul<&'a Iso2> for Iso2 {
    type Output = Iso2;
    #[inline]
    fn mul(self, rhs: &'a Iso2) -> Iso2 {
        &self * rhs
    }
}

impl<'a> Mul<Iso2> for &'a Iso2 {
    type Output = Iso2;
    #[inline]
    fn mul(self, rhs: Iso2) -> Iso2 {
        self * &rhs
    }
}

impl<'a, 'b> Mul<&'b Iso2> for &'a Iso2 {
    type Output = Iso2;
    #[inline]
    fn mul(self, rhs: &'b Iso2) -> Iso2 {
        Iso2 {
            rotation: self.rotation * rhs.rotation,
            translation: self.translation + self.rotation * rhs.translation,
        }
    }
}

impl IntoInner for Mat2 {}
impl IntoInner for Vec2 {}

impl MatrixOps<2> for Mat2 {
    type Vector = Vec2;
    type SymmetricEigen = SymmetricEigen2;

    #[inline(always)]
    fn zeros() -> Self {
        Self::ZERO
    }

    #[inline]
    fn from_columns(columns: &[Vec2; 2]) -> Self {
        Mat2::from_cols(columns[0], columns[1])
    }

    #[inline(always)]
    fn column_mut(&mut self, col: usize) -> &mut Self::Vector {
        self.col_mut(col)
    }

    #[cfg(feature = "dim2")]
    #[inline(always)]
    fn new(m00: f32, m01: f32, m10: f32, m11: f32) -> Self {
        Mat2::from_cols_array(&[m00, m10, m01, m11])
    }

    #[cfg(feature = "dim3")]
    fn new(
        _m00: f32,
        _m01: f32,
        _m02: f32,
        _m10: f32,
        _m11: f32,
        _m12: f32,
        _m20: f32,
        _m21: f32,
        _m22: f32,
    ) -> Self {
        unreachable!()
    }

    #[inline]
    fn from_diagonal_element(elt: Real) -> Self {
        Self::from_diagonal(Self::Vector::splat(elt))
    }

    #[inline]
    fn try_inverse(&self) -> Option<Self> {
        (self.determinant() != 0.0).then_some(self.inverse())
    }

    #[inline]
    fn abs(&self) -> Self {
        Self::from_cols(self.x_axis.abs(), self.y_axis.abs())
    }

    #[inline]
    fn symmetric_eigen(&self) -> SymmetricEigen2 {
        let na_mat = na::Matrix2::from(*self);
        let eig = na_mat.symmetric_eigen();
        SymmetricEigen2 {
            eigenvalues: eig.eigenvalues.into(),
            eigenvectors: eig.eigenvectors.into(),
        }
    }

    #[inline]
    fn symmetric_eigenvalues(&self) -> Self::Vector {
        let na_mat = na::Matrix2::from(*self);
        na_mat.symmetric_eigenvalues().into()
    }

    #[inline]
    fn swap_columns(&mut self, i: usize, j: usize) {
        let mut cols = [self.x_axis, self.y_axis];
        cols.swap(i, j);
        *self = Self::from_cols(cols[0], cols[1]);
    }

    #[inline]
    fn column(&self, i: usize) -> Self::Vector {
        self.col(i)
    }
}

impl PointOps for Vec2 {
    type Coords = Self;

    #[inline(always)]
    fn into_vector(self) -> Self {
        self
    }

    #[inline(always)]
    fn as_vector_mut(&mut self) -> &mut Self {
        self
    }

    #[inline(always)]
    fn origin() -> Self {
        Self::ZERO
    }
}

impl WSign<Vec2> for Vec2 {
    #[inline(always)]
    fn copy_sign_to(self, to: Vec2) -> Vec2 {
        to.copysign(self)
    }
}

impl GlamVectorOps<2> for Vec2 {
    type Matrix = Mat2;

    #[inline(always)]
    fn x() -> Self {
        Self::X
    }

    #[inline(always)]
    fn y() -> Self {
        Self::Y
    }

    #[cfg(feature = "dim3")]
    #[inline(always)]
    fn z() -> Self {
        unreachable!()
    }

    #[inline(always)]
    fn ith(i: usize, val: Real) -> Self {
        let mut result = Self::ZERO;
        result[i] = val;
        result
    }

    #[inline(always)]
    fn from_row_slice(elts: &[Real]) -> Self {
        Self::new(elts[0], elts[1])
    }

    #[inline(always)]
    fn repeat(value: Real) -> Self {
        Self::splat(value)
    }

    #[inline]
    fn orthonormal_subspace_basis(vs: &[Self], mut f: impl FnMut(&Self) -> bool) {
        // This is extracted from nalgebra’s generic implementation of orthonormal_subspace_basis.
        if vs.is_empty() {
            let _ = f(&Vec2::X) && f(&Vec2::Y);
        } else if vs.len() == 1 {
            let v = &vs[0];
            let res = Vec2::new(-v[1], v[0]);

            let _ = f(&res.normalize());
        }

        // Otherwise, nothing.
    }

    #[inline]
    fn iter(&self) -> std::array::IntoIter<&Real, 2> {
        [&self.x, &self.y].into_iter()
    }

    #[inline(always)]
    fn as_slice(&self) -> &[Real] {
        &self.as_ref()[..]
    }

    #[inline(always)]
    fn swap_rows(&mut self, i: usize, j: usize) {
        self.as_mut().swap(i, j);
    }

    #[inline(always)]
    fn map(&self, mut f: impl FnMut(Real) -> Real) -> Self {
        Self::new(f(self.x), f(self.y))
    }

    #[inline(always)]
    fn is_zero(self) -> bool {
        self == Self::ZERO
    }

    #[inline(always)]
    fn norm_squared(&self) -> Real {
        self.length_squared()
    }

    #[inline(always)]
    fn norm(&self) -> Real {
        self.length()
    }

    #[inline(always)]
    fn inf(&self, rhs: &Self) -> Self {
        self.min(*rhs)
    }

    #[inline(always)]
    fn sup(&self, rhs: &Self) -> Self {
        self.max(*rhs)
    }

    #[inline(always)]
    fn component_div(&self, rhs: &Self) -> Self {
        *self / *rhs
    }

    #[inline(always)]
    fn component_mul(&self, rhs: &Self) -> Self {
        *self * *rhs
    }

    #[inline(always)]
    fn component_mul_assign(&mut self, rhs: &Self) {
        *self *= *rhs;
    }

    #[inline(always)]
    fn zeros() -> Self {
        Self::ZERO
    }

    #[inline]
    fn normalize_mut(&mut self) -> Real {
        let norm = self.length();
        *self /= norm;
        norm
    }

    #[inline]
    fn try_normalize_mut(&mut self, eps: Real) -> Option<Real> {
        let norm = self.length();
        if norm > eps {
            *self /= norm;
            Some(norm)
        } else {
            None
        }
    }

    #[inline]
    fn imin(&self) -> usize {
        if self.x <= self.y {
            0
        } else {
            1
        }
    }

    #[inline]
    fn amin(&self) -> Real {
        self.x.abs().min(self.y.abs())
    }

    #[inline]
    fn amax(&self) -> Real {
        self.x.abs().max(self.y.abs())
    }

    #[inline]
    fn iamin(&self) -> usize {
        self.abs().imin()
    }

    #[inline]
    fn iamax(&self) -> usize {
        if self.x.abs() >= self.y.abs() {
            0
        } else {
            1
        }
    }

    #[rustfmt::skip]
    #[inline]
    #[cfg(feature = "dim2")]
    fn outer_product(&self, rhs: &Self) -> Mat2 {
        Mat2::new(
            self.x * rhs.x, self.x * rhs.y, self.y * rhs.x, self.y * rhs.y
        )
    }


    #[rustfmt::skip]
    #[inline]
    #[cfg(feature = "dim3")]
    fn outer_product(&self, _rhs: &Self) -> Mat2 {
        unreachable!()
    }

    #[inline(always)]
    fn angle(&self, rhs: &Self) -> Real {
        self.angle_between(*rhs)
    }

    #[inline(always)]
    fn zip_apply<F>(&mut self, b: &Self, mut f: F)
    where
        F: FnMut(&mut Real, Real),
    {
        f(&mut self.x, b.x);
        f(&mut self.y, b.y);
    }

    /*
     * Unit vector.
     */
    #[inline]
    fn try_new_and_get(v: Self, eps: Real) -> Option<(Self, Real)> {
        let sqnorm = v.norm_squared();
        if sqnorm > eps * eps {
            let norm = sqnorm.sqrt();
            Some((v / norm, norm))
        } else {
            None
        }
    }
}

impl WBasis for Vec2 {
    type Basis = [Vec2; 1];
    #[inline]
    fn orthonormal_basis(self) -> [Vec2; 1] {
        [Vec2::new(-self.y, self.x)]
    }
}

impl RotationOps for Mat2 {
    type Vector = Vec2;
    type UnitVector = Vec2;
    type AngVector = Real;
    type RotationMatrix = Mat2;
    type RotationBetween = Self;

    #[inline(always)]
    fn identity() -> Self {
        Self::IDENTITY
    }

    #[inline(always)]
    fn from_scaled_axis(axis_angle: Self::AngVector) -> Self {
        Mat2::from_angle(axis_angle)
    }

    #[inline(always)]
    fn from_rotation_matrix(mat: &Self::RotationMatrix) -> Self {
        *mat
    }

    #[inline(always)]
    fn to_rotation_matrix(&self) -> Self::RotationMatrix {
        *self
    }

    #[inline(always)]
    fn renormalize(&mut self) {
        let norm = self.x_axis.norm();
        *self *= 1.0 / norm;
    }

    #[inline(always)]
    fn renormalize_fast(&mut self) {
        self.renormalize();
    }

    #[inline(always)]
    fn new_eps(axis_angle: Self::AngVector, val: Real) -> Self {
        if axis_angle.abs() <= val {
            Self::identity()
        } else {
            Self::from_angle(axis_angle)
        }
    }

    #[inline(always)]
    fn rotation_between(a: &Self::Vector, b: &Self::Vector) -> Self {
        Mat2::from_angle(a.angle_between(*b))
    }

    #[inline(always)]
    fn scaled_rotation_between_axis(na: &Self::Vector, nb: &Self::Vector, s: Real) -> Self {
        let sang = na.perp_dot(*nb);
        let cang = na.dot(*nb);
        Self::from_angle(sang.simd_atan2(cang) * s)
    }

    #[inline(always)]
    fn scaled_axis(&self) -> Self::AngVector {
        unreachable!()
    }

    #[inline(always)]
    fn axis_angle(&self) -> Option<(Self::UnitVector, Real)> {
        unreachable!()
    }

    #[inline(always)]
    fn imag(&self) -> Self::AngVector {
        self.x_axis.y
    }
}

impl WCross<Vec2> for Vec2 {
    type Result = Real;
    #[inline]
    fn gcross(&self, rhs: Vec2) -> Self::Result {
        self.perp_dot(rhs)
    }
}

impl WCross<Vec2> for Real {
    type Result = Vec2;

    fn gcross(&self, rhs: Vec2) -> Self::Result {
        Vec2::new(-rhs.y * *self, rhs.x * *self)
    }
}

impl VectorOps for Vec2 {
    #[inline]
    fn try_normalize_eps(mut self, eps: Real) -> Option<Self> {
        let result = self.try_normalize_mut(eps);
        result.map(|_| self)
    }
}

impl VectorU32Ops for UVec2 {
    type Vector = Vec2;

    #[inline]
    fn zeros() -> Self {
        UVec2::ZERO
    }

    #[inline]
    fn repeat(val: u32) -> Self {
        UVec2::splat(val)
    }

    #[inline]
    fn inf(&self, rhs: &Self) -> Self {
        self.min(*rhs)
    }

    #[inline]
    fn sup(&self, rhs: &Self) -> Self {
        self.max(*rhs)
    }

    #[inline(always)]
    fn map(&self, mut f: impl FnMut(u32) -> Real) -> Vec2 {
        Vec2::new(f(self.x), f(self.y))
    }

    #[inline(always)]
    fn apply(&mut self, mut f: impl FnMut(&mut u32)) {
        f(&mut self.x);
        f(&mut self.y);
    }
}

// impl From<Iso2> for Isometry2<Real> {
//     #[inline]
//     fn from(value: Iso2) -> Self {
//         (value.translation, value.rotation).into()
//     }
// }

impl Mul<Vec2> for SdpMatrix2<Real> {
    type Output = Vec2;

    #[inline]
    fn mul(self, rhs: Vec2) -> Self::Output {
        let x = self.m11 * rhs.x + self.m12 * rhs.y;
        let y = self.m12 * rhs.x + self.m22 * rhs.y;
        Vec2::new(x, y)
    }
}

impl Mul<Mat2> for SdpMatrix2<Real> {
    type Output = Mat2;

    #[inline]
    fn mul(self, rhs: Mat2) -> Self::Output {
        let x0 = self.m11 * rhs.x_axis.x + self.m12 * rhs.x_axis.y;
        let y0 = self.m12 * rhs.x_axis.x + self.m22 * rhs.x_axis.y;

        let x1 = self.m11 * rhs.y_axis.x + self.m12 * rhs.y_axis.y;
        let y1 = self.m12 * rhs.y_axis.x + self.m22 * rhs.y_axis.y;

        Mat2::from_cols_array(&[x0, y0, x1, y1])
    }
}
