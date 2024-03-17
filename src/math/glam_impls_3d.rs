pub use super::{
    GlamVectorOps, MatrixOps, PointOps, Real, RotationOps, SymmetricEigen3, VectorOps, VectorU32Ops,
};
use crate::math::IntoInner;
use crate::utils::{SdpMatrix3, WBasis, WSign};
use core::ops::Mul;
use glam::{Mat3, Quat, UVec3, Vec3};
use na::Isometry3;

/// A 3D isometry.
#[derive(Copy, Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde-serialize", derive(Serialize, Deserialize))]
pub struct Iso3 {
    pub rotation: Quat,
    pub translation: Vec3,
}

impl Default for Iso3 {
    #[inline]
    fn default() -> Self {
        Self::identity()
    }
}

impl Iso3 {
    #[inline]
    pub fn new(translation: Vec3, axis_angle: Vec3) -> Self {
        Self {
            rotation: Quat::from_scaled_axis(axis_angle),
            translation,
        }
    }

    #[inline]
    pub fn from_parts(translation: Vec3, rotation: Quat) -> Self {
        Self {
            rotation,
            translation,
        }
    }

    #[inline]
    pub fn identity() -> Self {
        Self {
            rotation: Quat::IDENTITY,
            translation: Vec3::ZERO,
        }
    }

    #[inline]
    pub fn rotation(axis_angle: Vec3) -> Self {
        Self {
            rotation: Quat::from_scaled_axis(axis_angle),
            translation: Vec3::ZERO,
        }
    }

    #[inline]
    pub fn translation(x: Real, y: Real, z: Real) -> Self {
        Self {
            rotation: Quat::IDENTITY,
            translation: Vec3::new(x, y, z),
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
        Self {
            rotation: self.rotation.slerp(rhs.rotation, t),
            translation: self.translation.lerp(rhs.translation, t),
        }
    }
}

impl From<Vec3> for Iso3 {
    #[inline]
    fn from(translation: Vec3) -> Self {
        Self {
            rotation: Quat::IDENTITY,
            translation,
        }
    }
}

impl Mul<Iso3> for Quat {
    type Output = Iso3;
    #[inline]
    fn mul(self, rhs: Iso3) -> Iso3 {
        Iso3 {
            rotation: self * rhs.rotation,
            translation: self * rhs.translation,
        }
    }
}

impl Mul<Iso3> for Iso3 {
    type Output = Iso3;
    #[inline]
    fn mul(self, rhs: Iso3) -> Iso3 {
        &self * &rhs
    }
}

impl<'a> Mul<&'a Iso3> for Iso3 {
    type Output = Iso3;
    #[inline]
    fn mul(self, rhs: &'a Iso3) -> Iso3 {
        &self * rhs
    }
}

impl<'a> Mul<Iso3> for &'a Iso3 {
    type Output = Iso3;
    #[inline]
    fn mul(self, rhs: Iso3) -> Iso3 {
        self * &rhs
    }
}

impl<'a, 'b> Mul<&'b Iso3> for &'a Iso3 {
    type Output = Iso3;
    #[inline]
    fn mul(self, rhs: &'b Iso3) -> Iso3 {
        Iso3 {
            rotation: self.rotation * rhs.rotation,
            translation: self.translation + self.rotation * rhs.translation,
        }
    }
}

impl MatrixOps<3> for Mat3 {
    type Vector = Vec3;
    type SymmetricEigen = SymmetricEigen3;

    #[inline(always)]
    fn zeros() -> Self {
        Self::ZERO
    }

    #[inline]
    fn from_columns(columns: &[Vec3; 3]) -> Self {
        Self {
            x_axis: columns[0],
            y_axis: columns[1],
            z_axis: columns[2],
        }
    }

    #[inline(always)]
    fn column_mut(&mut self, col: usize) -> &mut Vec3 {
        self.col_mut(col)
    }

    #[inline(always)]
    #[cfg(feature = "dim2")]
    fn new(m00: f32, m01: f32, m10: f32, m11: f32) -> Self {
        Mat3::from_cols_array(&[m00, m10, 0.0, m01, m11, 0.0, 0.0, 0.0, 1.0])
    }

    #[inline(always)]
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
    ) -> Self {
        Mat3::from_cols_array(&[m00, m10, m20, m01, m11, m21, m02, m12, m22])
    }

    #[inline]
    fn from_diagonal_element(elt: Real) -> Self {
        Self::from_diagonal(Vec3::splat(elt))
    }

    #[inline]
    fn try_inverse(&self) -> Option<Self> {
        (self.determinant() != 0.0).then_some(self.inverse())
    }

    #[inline]
    fn abs(&self) -> Self {
        Self::from_cols(self.x_axis.abs(), self.y_axis.abs(), self.z_axis.abs())
    }

    #[inline]
    fn symmetric_eigen(&self) -> Self::SymmetricEigen {
        let na_mat = na::Matrix3::from(*self);
        let eig = na_mat.symmetric_eigen();
        SymmetricEigen3 {
            eigenvalues: eig.eigenvalues.into(),
            eigenvectors: eig.eigenvectors.into(),
        }
    }

    #[inline]
    fn symmetric_eigenvalues(&self) -> Vec3 {
        let na_mat = na::Matrix3::from(*self);
        na_mat.symmetric_eigenvalues().into()
    }

    #[inline]
    fn swap_columns(&mut self, i: usize, j: usize) {
        let mut cols = [self.x_axis, self.y_axis, self.z_axis];
        cols.swap(i, j);
        *self = Self::from_cols(cols[0], cols[1], cols[2]);
    }

    #[inline]
    fn column(&self, i: usize) -> Vec3 {
        self.col(i)
    }
}

impl PointOps for Vec3 {
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

impl WSign<Vec3> for Vec3 {
    #[inline(always)]
    fn copy_sign_to(self, to: Vec3) -> Vec3 {
        to.copysign(self)
    }
}

impl GlamVectorOps<3> for Vec3 {
    type Matrix = Mat3;

    #[inline(always)]
    fn x() -> Self {
        Self::X
    }

    #[inline(always)]
    fn y() -> Self {
        Self::Y
    }

    #[inline(always)]
    #[cfg(feature = "dim3")]
    fn z() -> Self {
        Self::Z
    }

    #[inline(always)]
    fn ith(i: usize, val: Real) -> Self {
        let mut result = Self::ZERO;
        result[i] = val;
        result
    }

    #[inline(always)]
    fn from_row_slice(elts: &[Real]) -> Self {
        Self::new(elts[0], elts[1], elts[2])
    }

    #[inline(always)]
    fn repeat(value: Real) -> Self {
        Self::splat(value)
    }

    #[inline(always)]
    fn swap_rows(&mut self, i: usize, j: usize) {
        self.as_mut().swap(i, j);
    }

    #[inline]
    fn orthonormal_subspace_basis(vs: &[Self], mut f: impl FnMut(&Self) -> bool) {
        // This is extracted from nalgebra’s generic implementation of orthonormal_subspace_basis.
        if vs.is_empty() {
            let _ = f(&Vec3::X) && f(&Vec3::Y) && f(&Vec3::Z);
        } else if vs.len() == 1 {
            let v = &vs[0];
            let mut a;

            if v[0].abs() > v[1].abs() {
                a = Vec3::new(v[2], 0.0, -v[0]);
            } else {
                a = Vec3::new(0.0, -v[2], v[1]);
            };

            let _ = a.normalize_mut();

            if f(&a.cross(*v)) {
                let _ = f(&a);
            }
        } else if vs.len() == 2 {
            let _ = f(&vs[0].cross(vs[1]).normalize());
        }
    }

    #[inline]
    fn iter(&self) -> std::array::IntoIter<&Real, 3> {
        [&self.x, &self.y, &self.z].into_iter()
    }

    #[inline(always)]
    fn map(&self, mut f: impl FnMut(Real) -> Real) -> Self {
        Vec3::new(f(self.x), f(self.y), f(self.z))
    }

    #[inline(always)]
    fn as_slice(&self) -> &[Real] {
        &self.as_ref()[..]
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
            if self.x <= self.z {
                0
            } else {
                2
            }
        } else if self.y <= self.z {
            1
        } else {
            2
        }
    }

    #[inline]
    fn amin(&self) -> Real {
        self.x.abs().min(self.y.abs()).min(self.z.abs())
    }

    #[inline]
    fn amax(&self) -> Real {
        self.x.abs().max(self.y.abs()).max(self.z.abs())
    }

    #[inline]
    fn iamin(&self) -> usize {
        self.abs().imin()
    }

    #[inline]
    fn iamax(&self) -> usize {
        let abs = self.abs();

        if abs.x >= abs.y {
            if abs.x >= abs.z {
                0
            } else {
                2
            }
        } else if abs.y >= abs.z {
            1
        } else {
            2
        }
    }

    #[rustfmt::skip]
    #[inline]
    fn outer_product(&self, rhs: &Self) -> Self::Matrix {
        Mat3::from_cols(
            *self * rhs.x, *self * rhs.y, *self * rhs.z
        )
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
        f(&mut self.z, b.z);
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

impl crate::utils::WCross<Vec3> for Vec3 {
    type Result = Vec3;
    #[inline]
    fn gcross(&self, rhs: Vec3) -> Self::Result {
        self.cross(rhs)
    }
}

impl WBasis for Vec3 {
    type Basis = [Vec3; 2];

    #[inline]
    fn orthonormal_basis(self) -> [Vec3; 2] {
        let sign = self.z.copy_sign_to(1.0);
        let a = -1.0 / (sign + self.z);
        let b = self.x * self.y * a;

        [
            Vec3::new(1.0 + sign * self.x * self.x * a, sign * b, -sign * self.x),
            Vec3::new(b, sign + self.y * self.y * a, -self.y),
        ]
    }
}

impl RotationOps for Quat {
    type Vector = Vec3;
    type UnitVector = Vec3;
    type AngVector = Vec3;
    type RotationMatrix = Mat3;
    type RotationBetween = Option<Self>;

    #[inline(always)]
    fn identity() -> Self {
        Self::IDENTITY
    }

    #[inline(always)]
    fn from_scaled_axis(axis_angle: Vec3) -> Self {
        Quat::from_scaled_axis(axis_angle)
    }

    #[inline(always)]
    fn from_rotation_matrix(mat: &Self::RotationMatrix) -> Self {
        Self::from_mat3(&mat)
    }

    #[inline(always)]
    fn to_rotation_matrix(&self) -> Self::RotationMatrix {
        Self::RotationMatrix::from_quat(*self)
    }

    #[inline(always)]
    fn renormalize(&mut self) {
        *self = self.normalize();
    }

    #[inline(always)]
    fn renormalize_fast(&mut self) {
        *self = self.normalize();
    }

    #[inline(always)]
    fn new_eps(axis_angle: Vec3, val: Real) -> Self {
        if axis_angle.norm() <= val {
            Self::identity()
        } else {
            Self::from_scaled_axis(axis_angle)
        }
    }

    #[inline(always)]
    fn rotation_between(a: &Vec3, b: &Vec3) -> Option<Self> {
        Some(Quat::from_rotation_arc(a.normalize(), b.normalize()))
    }

    #[inline(always)]
    fn scaled_rotation_between_axis(na: &Vec3, nb: &Vec3, s: Real) -> Option<Self> {
        use crate::approx::AbsDiffEq;

        // Copied and adapted from nalgebra.
        let c = na.cross(*nb);

        if let Some(axis) = c.try_normalize_eps(Real::default_epsilon()) {
            let cos = na.dot(*nb);

            // The cosinus may be out of [-1, 1] because of inaccuracies.
            if cos <= -1.0 {
                None
            } else if cos >= 1.0 {
                Some(Self::identity())
            } else {
                Some(Self::from_axis_angle(axis, cos.acos() * s))
            }
        } else if na.dot(*nb) < 0.0 {
            // PI
            //
            // The rotation axis is undefined but the angle not zero. This is not a
            // simple rotation.
            None
        } else {
            // Zero
            Some(Self::identity())
        }
    }

    #[inline(always)]
    fn scaled_axis(&self) -> Self::AngVector {
        self.to_scaled_axis()
    }

    #[inline(always)]
    fn axis_angle(&self) -> Option<(Self::UnitVector, Real)> {
        Some(self.to_axis_angle())
    }

    #[inline(always)]
    fn imag(&self) -> Self::AngVector {
        self.xyz()
    }
}

impl VectorOps for Vec3 {
    #[inline]
    fn try_normalize_eps(mut self, eps: Real) -> Option<Self> {
        let result = self.try_normalize_mut(eps);
        result.map(|_| self)
    }
}

impl VectorU32Ops for UVec3 {
    type Vector = Vec3;

    #[inline]
    fn zeros() -> Self {
        UVec3::ZERO
    }

    #[inline]
    fn repeat(val: u32) -> Self {
        UVec3::splat(val)
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
    fn map(&self, mut f: impl FnMut(u32) -> Real) -> Vec3 {
        Vec3::new(f(self.x), f(self.y), f(self.z))
    }

    #[inline(always)]
    fn apply(&mut self, mut f: impl FnMut(&mut u32)) {
        f(&mut self.x);
        f(&mut self.y);
        f(&mut self.z);
    }
}

impl From<Iso3> for Isometry3<Real> {
    #[inline]
    fn from(value: Iso3) -> Self {
        (value.translation, value.rotation).into()
    }
}

impl Mul<Vec3> for SdpMatrix3<Real> {
    type Output = Vec3;

    #[inline]
    fn mul(self, rhs: Vec3) -> Self::Output {
        let x = self.m11 * rhs.x + self.m12 * rhs.y + self.m13 * rhs.z;
        let y = self.m12 * rhs.x + self.m22 * rhs.y + self.m23 * rhs.z;
        let z = self.m13 * rhs.x + self.m23 * rhs.y + self.m33 * rhs.z;
        Vec3::new(x, y, z)
    }
}

#[cfg(feature = "dim3")]
impl Mul<Mat3> for SdpMatrix3<Real> {
    type Output = Mat3;

    #[inline]
    fn mul(self, rhs: Mat3) -> Self::Output {
        let x0 = self.m11 * rhs.x_axis.x + self.m12 * rhs.x_axis.y + self.m13 * rhs.x_axis.z;
        let y0 = self.m12 * rhs.x_axis.x + self.m22 * rhs.x_axis.y + self.m23 * rhs.x_axis.z;
        let z0 = self.m13 * rhs.x_axis.x + self.m23 * rhs.x_axis.y + self.m33 * rhs.x_axis.z;

        let x1 = self.m11 * rhs.y_axis.x + self.m12 * rhs.y_axis.y + self.m13 * rhs.y_axis.z;
        let y1 = self.m12 * rhs.y_axis.x + self.m22 * rhs.y_axis.y + self.m23 * rhs.y_axis.z;
        let z1 = self.m13 * rhs.y_axis.x + self.m23 * rhs.y_axis.y + self.m33 * rhs.y_axis.z;

        let x2 = self.m11 * rhs.z_axis.x + self.m12 * rhs.z_axis.y + self.m13 * rhs.z_axis.z;
        let y2 = self.m12 * rhs.z_axis.x + self.m22 * rhs.z_axis.y + self.m23 * rhs.z_axis.z;
        let z2 = self.m13 * rhs.z_axis.x + self.m23 * rhs.z_axis.y + self.m33 * rhs.z_axis.z;

        Mat3::new(x0, x1, x2, y0, y1, y2, z0, z1, z2)
    }
}

impl IntoInner for Mat3 {}
impl IntoInner for Vec3 {}
impl IntoInner for Quat {}
