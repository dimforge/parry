use crate::math::Real;
use na::{Matrix2, Matrix3, Matrix3x2, SimdRealField, Vector2, Vector3};
use std::ops::{Add, Mul};

/// A 2x2 symmetric-definite-positive matrix.
#[derive(Copy, Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde-serialize", derive(Serialize, Deserialize))]
pub struct SdpMatrix2<N> {
    /// The component at the first row and first column of this matrix.
    pub m11: N,
    /// The component at the first row and second column of this matrix.
    pub m12: N,
    /// The component at the second row and second column of this matrix.
    pub m22: N,
}

impl<N: SimdRealField + Copy> SdpMatrix2<N> {
    /// A new SDP 2x2 matrix with the given components.
    ///
    /// Because the matrix is symmetric, only the lower off-diagonal component is required.
    pub fn new(m11: N, m12: N, m22: N) -> Self {
        Self { m11, m12, m22 }
    }

    /// Build an `SdpMatrix2` structure from a plain matrix, assuming it is SDP.
    ///
    /// No check is performed to ensure `mat` is actually SDP.
    pub fn from_sdp_matrix(mat: Matrix2<N>) -> Self {
        Self {
            m11: mat.m11,
            m12: mat.m12,
            m22: mat.m22,
        }
    }

    /// Create a new SDP matrix filled with zeros.
    pub fn zero() -> Self {
        Self {
            m11: N::zero(),
            m12: N::zero(),
            m22: N::zero(),
        }
    }

    /// Create a new SDP matrix with its diagonal filled with `val`, and its off-diagonal elements set to zero.
    pub fn diagonal(val: N) -> Self {
        Self {
            m11: val,
            m12: N::zero(),
            m22: val,
        }
    }

    /// Adds `val` to the diagonal components of `self`.
    pub fn add_diagonal(&mut self, elt: N) -> Self {
        Self {
            m11: self.m11 + elt,
            m12: self.m12,
            m22: self.m22 + elt,
        }
    }

    /// Compute the inverse of this SDP matrix without performing any inversibility check.
    pub fn inverse_unchecked(&self) -> Self {
        let determinant = self.m11 * self.m22 - self.m12 * self.m12;
        let m11 = self.m22 / determinant;
        let m12 = -self.m12 / determinant;
        let m22 = self.m11 / determinant;

        Self { m11, m12, m22 }
    }

    /// Convert this SDP matrix to a regular matrix representation.
    pub fn into_matrix(self) -> Matrix2<N> {
        Matrix2::new(self.m11, self.m12, self.m12, self.m22)
    }
}

impl<N: SimdRealField + Copy> Add<SdpMatrix2<N>> for SdpMatrix2<N> {
    type Output = Self;

    fn add(self, rhs: SdpMatrix2<N>) -> Self {
        Self::new(self.m11 + rhs.m11, self.m12 + rhs.m12, self.m22 + rhs.m22)
    }
}

impl<N: SimdRealField + Copy> Mul<Vector2<N>> for SdpMatrix2<N> {
    type Output = Vector2<N>;

    fn mul(self, rhs: Vector2<N>) -> Self::Output {
        Vector2::new(
            self.m11 * rhs.x + self.m12 * rhs.y,
            self.m12 * rhs.x + self.m22 * rhs.y,
        )
    }
}

impl Mul<Real> for SdpMatrix2<Real> {
    type Output = SdpMatrix2<Real>;

    fn mul(self, rhs: Real) -> Self::Output {
        SdpMatrix2::new(self.m11 * rhs, self.m12 * rhs, self.m22 * rhs)
    }
}

/// A 3x3 symmetric-definite-positive matrix.
#[derive(Copy, Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde-serialize", derive(Serialize, Deserialize))]
pub struct SdpMatrix3<N> {
    /// The component at the first row and first column of this matrix.
    pub m11: N,
    /// The component at the first row and second column of this matrix.
    pub m12: N,
    /// The component at the first row and third column of this matrix.
    pub m13: N,
    /// The component at the second row and second column of this matrix.
    pub m22: N,
    /// The component at the second row and third column of this matrix.
    pub m23: N,
    /// The component at the third row and third column of this matrix.
    pub m33: N,
}

impl<N: SimdRealField + Copy> SdpMatrix3<N> {
    /// A new SDP 3x3 matrix with the given components.
    ///
    /// Because the matrix is symmetric, only the lower off-diagonal components is required.
    pub fn new(m11: N, m12: N, m13: N, m22: N, m23: N, m33: N) -> Self {
        Self {
            m11,
            m12,
            m13,
            m22,
            m23,
            m33,
        }
    }

    /// Build an `SdpMatrix3` structure from a plain matrix, assuming it is SDP.
    ///
    /// No check is performed to ensure `mat` is actually SDP.
    pub fn from_sdp_matrix(mat: Matrix3<N>) -> Self {
        Self {
            m11: mat.m11,
            m12: mat.m12,
            m13: mat.m13,
            m22: mat.m22,
            m23: mat.m23,
            m33: mat.m33,
        }
    }

    /// Create a new SDP matrix filled with zeros.
    pub fn zero() -> Self {
        Self {
            m11: N::zero(),
            m12: N::zero(),
            m13: N::zero(),
            m22: N::zero(),
            m23: N::zero(),
            m33: N::zero(),
        }
    }

    /// Create a new SDP matrix with its diagonal filled with `val`, and its off-diagonal elements set to zero.
    pub fn diagonal(val: N) -> Self {
        Self {
            m11: val,
            m12: N::zero(),
            m13: N::zero(),
            m22: val,
            m23: N::zero(),
            m33: val,
        }
    }

    /// Are all components of this matrix equal to zero?
    pub fn is_zero(&self) -> bool {
        self.m11.is_zero()
            && self.m12.is_zero()
            && self.m13.is_zero()
            && self.m22.is_zero()
            && self.m23.is_zero()
            && self.m33.is_zero()
    }

    /// Compute the inverse of this SDP matrix without performing any inversibility check.
    pub fn inverse_unchecked(&self) -> Self {
        let minor_m12_m23 = self.m22 * self.m33 - self.m23 * self.m23;
        let minor_m11_m23 = self.m12 * self.m33 - self.m13 * self.m23;
        let minor_m11_m22 = self.m12 * self.m23 - self.m13 * self.m22;

        let determinant =
            self.m11 * minor_m12_m23 - self.m12 * minor_m11_m23 + self.m13 * minor_m11_m22;
        let inv_det = N::one() / determinant;

        SdpMatrix3 {
            m11: minor_m12_m23 * inv_det,
            m12: -minor_m11_m23 * inv_det,
            m13: minor_m11_m22 * inv_det,
            m22: (self.m11 * self.m33 - self.m13 * self.m13) * inv_det,
            m23: (self.m13 * self.m12 - self.m23 * self.m11) * inv_det,
            m33: (self.m11 * self.m22 - self.m12 * self.m12) * inv_det,
        }
    }

    /// Compute the quadratic form `m.transpose() * self * m`.
    pub fn quadform3x2(&self, m: &Matrix3x2<N>) -> SdpMatrix2<N> {
        let x0 = self.m11 * m.m11 + self.m12 * m.m21 + self.m13 * m.m31;
        let y0 = self.m12 * m.m11 + self.m22 * m.m21 + self.m23 * m.m31;
        let z0 = self.m13 * m.m11 + self.m23 * m.m21 + self.m33 * m.m31;

        let x1 = self.m11 * m.m12 + self.m12 * m.m22 + self.m13 * m.m32;
        let y1 = self.m12 * m.m12 + self.m22 * m.m22 + self.m23 * m.m32;
        let z1 = self.m13 * m.m12 + self.m23 * m.m22 + self.m33 * m.m32;

        let m11 = m.m11 * x0 + m.m21 * y0 + m.m31 * z0;
        let m12 = m.m11 * x1 + m.m21 * y1 + m.m31 * z1;
        let m22 = m.m12 * x1 + m.m22 * y1 + m.m32 * z1;

        SdpMatrix2 { m11, m12, m22 }
    }

    /// Compute the quadratic form `m.transpose() * self * m`.
    pub fn quadform(&self, m: &Matrix3<N>) -> Self {
        let x0 = self.m11 * m.m11 + self.m12 * m.m21 + self.m13 * m.m31;
        let y0 = self.m12 * m.m11 + self.m22 * m.m21 + self.m23 * m.m31;
        let z0 = self.m13 * m.m11 + self.m23 * m.m21 + self.m33 * m.m31;

        let x1 = self.m11 * m.m12 + self.m12 * m.m22 + self.m13 * m.m32;
        let y1 = self.m12 * m.m12 + self.m22 * m.m22 + self.m23 * m.m32;
        let z1 = self.m13 * m.m12 + self.m23 * m.m22 + self.m33 * m.m32;

        let x2 = self.m11 * m.m13 + self.m12 * m.m23 + self.m13 * m.m33;
        let y2 = self.m12 * m.m13 + self.m22 * m.m23 + self.m23 * m.m33;
        let z2 = self.m13 * m.m13 + self.m23 * m.m23 + self.m33 * m.m33;

        let m11 = m.m11 * x0 + m.m21 * y0 + m.m31 * z0;
        let m12 = m.m11 * x1 + m.m21 * y1 + m.m31 * z1;
        let m13 = m.m11 * x2 + m.m21 * y2 + m.m31 * z2;

        let m22 = m.m12 * x1 + m.m22 * y1 + m.m32 * z1;
        let m23 = m.m12 * x2 + m.m22 * y2 + m.m32 * z2;
        let m33 = m.m13 * x2 + m.m23 * y2 + m.m33 * z2;

        Self {
            m11,
            m12,
            m13,
            m22,
            m23,
            m33,
        }
    }

    /// Adds `elt` to the diagonal components of `self`.
    pub fn add_diagonal(&self, elt: N) -> Self {
        Self {
            m11: self.m11 + elt,
            m12: self.m12,
            m13: self.m13,
            m22: self.m22 + elt,
            m23: self.m23,
            m33: self.m33 + elt,
        }
    }
}

impl<N: Add<N>> Add<SdpMatrix3<N>> for SdpMatrix3<N> {
    type Output = SdpMatrix3<N::Output>;

    fn add(self, rhs: SdpMatrix3<N>) -> Self::Output {
        SdpMatrix3 {
            m11: self.m11 + rhs.m11,
            m12: self.m12 + rhs.m12,
            m13: self.m13 + rhs.m13,
            m22: self.m22 + rhs.m22,
            m23: self.m23 + rhs.m23,
            m33: self.m33 + rhs.m33,
        }
    }
}

impl Mul<Real> for SdpMatrix3<Real> {
    type Output = SdpMatrix3<Real>;

    fn mul(self, rhs: Real) -> Self::Output {
        SdpMatrix3 {
            m11: self.m11 * rhs,
            m12: self.m12 * rhs,
            m13: self.m13 * rhs,
            m22: self.m22 * rhs,
            m23: self.m23 * rhs,
            m33: self.m33 * rhs,
        }
    }
}

impl<N: SimdRealField + Copy> Mul<Vector3<N>> for SdpMatrix3<N> {
    type Output = Vector3<N>;

    fn mul(self, rhs: Vector3<N>) -> Self::Output {
        let x = self.m11 * rhs.x + self.m12 * rhs.y + self.m13 * rhs.z;
        let y = self.m12 * rhs.x + self.m22 * rhs.y + self.m23 * rhs.z;
        let z = self.m13 * rhs.x + self.m23 * rhs.y + self.m33 * rhs.z;
        Vector3::new(x, y, z)
    }
}

impl<N: SimdRealField + Copy> Mul<Matrix3<N>> for SdpMatrix3<N> {
    type Output = Matrix3<N>;

    fn mul(self, rhs: Matrix3<N>) -> Self::Output {
        let x0 = self.m11 * rhs.m11 + self.m12 * rhs.m21 + self.m13 * rhs.m31;
        let y0 = self.m12 * rhs.m11 + self.m22 * rhs.m21 + self.m23 * rhs.m31;
        let z0 = self.m13 * rhs.m11 + self.m23 * rhs.m21 + self.m33 * rhs.m31;

        let x1 = self.m11 * rhs.m12 + self.m12 * rhs.m22 + self.m13 * rhs.m32;
        let y1 = self.m12 * rhs.m12 + self.m22 * rhs.m22 + self.m23 * rhs.m32;
        let z1 = self.m13 * rhs.m12 + self.m23 * rhs.m22 + self.m33 * rhs.m32;

        let x2 = self.m11 * rhs.m13 + self.m12 * rhs.m23 + self.m13 * rhs.m33;
        let y2 = self.m12 * rhs.m13 + self.m22 * rhs.m23 + self.m23 * rhs.m33;
        let z2 = self.m13 * rhs.m13 + self.m23 * rhs.m23 + self.m33 * rhs.m33;

        Matrix3::new(x0, x1, x2, y0, y1, y2, z0, z1, z2)
    }
}

impl<N: SimdRealField + Copy> Mul<Matrix3x2<N>> for SdpMatrix3<N> {
    type Output = Matrix3x2<N>;

    fn mul(self, rhs: Matrix3x2<N>) -> Self::Output {
        let x0 = self.m11 * rhs.m11 + self.m12 * rhs.m21 + self.m13 * rhs.m31;
        let y0 = self.m12 * rhs.m11 + self.m22 * rhs.m21 + self.m23 * rhs.m31;
        let z0 = self.m13 * rhs.m11 + self.m23 * rhs.m21 + self.m33 * rhs.m31;

        let x1 = self.m11 * rhs.m12 + self.m12 * rhs.m22 + self.m13 * rhs.m32;
        let y1 = self.m12 * rhs.m12 + self.m22 * rhs.m22 + self.m23 * rhs.m32;
        let z1 = self.m13 * rhs.m12 + self.m23 * rhs.m22 + self.m33 * rhs.m32;

        Matrix3x2::new(x0, x1, y0, y1, z0, z1)
    }
}

impl<T> From<[SdpMatrix3<Real>; 4]> for SdpMatrix3<T>
where
    T: From<[Real; 4]>,
{
    fn from(data: [SdpMatrix3<Real>; 4]) -> Self {
        SdpMatrix3 {
            m11: T::from([data[0].m11, data[1].m11, data[2].m11, data[3].m11]),
            m12: T::from([data[0].m12, data[1].m12, data[2].m12, data[3].m12]),
            m13: T::from([data[0].m13, data[1].m13, data[2].m13, data[3].m13]),
            m22: T::from([data[0].m22, data[1].m22, data[2].m22, data[3].m22]),
            m23: T::from([data[0].m23, data[1].m23, data[2].m23, data[3].m23]),
            m33: T::from([data[0].m33, data[1].m33, data[2].m33, data[3].m33]),
        }
    }
}

#[cfg(feature = "simd-nightly")]
impl From<[SdpMatrix3<Real>; 8]> for SdpMatrix3<simba::simd::f32x8> {
    fn from(data: [SdpMatrix3<Real>; 8]) -> Self {
        SdpMatrix3 {
            m11: simba::simd::f32x8::from([
                data[0].m11,
                data[1].m11,
                data[2].m11,
                data[3].m11,
                data[4].m11,
                data[5].m11,
                data[6].m11,
                data[7].m11,
            ]),
            m12: simba::simd::f32x8::from([
                data[0].m12,
                data[1].m12,
                data[2].m12,
                data[3].m12,
                data[4].m12,
                data[5].m12,
                data[6].m12,
                data[7].m12,
            ]),
            m13: simba::simd::f32x8::from([
                data[0].m13,
                data[1].m13,
                data[2].m13,
                data[3].m13,
                data[4].m13,
                data[5].m13,
                data[6].m13,
                data[7].m13,
            ]),
            m22: simba::simd::f32x8::from([
                data[0].m22,
                data[1].m22,
                data[2].m22,
                data[3].m22,
                data[4].m22,
                data[5].m22,
                data[6].m22,
                data[7].m22,
            ]),
            m23: simba::simd::f32x8::from([
                data[0].m23,
                data[1].m23,
                data[2].m23,
                data[3].m23,
                data[4].m23,
                data[5].m23,
                data[6].m23,
                data[7].m23,
            ]),
            m33: simba::simd::f32x8::from([
                data[0].m33,
                data[1].m33,
                data[2].m33,
                data[3].m33,
                data[4].m33,
                data[5].m33,
                data[6].m33,
                data[7].m33,
            ]),
        }
    }
}

#[cfg(feature = "simd-nightly")]
impl From<[SdpMatrix3<Real>; 16]> for SdpMatrix3<simba::simd::f32x16> {
    fn from(data: [SdpMatrix3<Real>; 16]) -> Self {
        SdpMatrix3 {
            m11: simba::simd::f32x16::from([
                data[0].m11,
                data[1].m11,
                data[2].m11,
                data[3].m11,
                data[4].m11,
                data[5].m11,
                data[6].m11,
                data[7].m11,
                data[8].m11,
                data[9].m11,
                data[10].m11,
                data[11].m11,
                data[12].m11,
                data[13].m11,
                data[14].m11,
                data[15].m11,
            ]),
            m12: simba::simd::f32x16::from([
                data[0].m12,
                data[1].m12,
                data[2].m12,
                data[3].m12,
                data[4].m12,
                data[5].m12,
                data[6].m12,
                data[7].m12,
                data[8].m12,
                data[9].m12,
                data[10].m12,
                data[11].m12,
                data[12].m12,
                data[13].m12,
                data[14].m12,
                data[15].m12,
            ]),
            m13: simba::simd::f32x16::from([
                data[0].m13,
                data[1].m13,
                data[2].m13,
                data[3].m13,
                data[4].m13,
                data[5].m13,
                data[6].m13,
                data[7].m13,
                data[8].m13,
                data[9].m13,
                data[10].m13,
                data[11].m13,
                data[12].m13,
                data[13].m13,
                data[14].m13,
                data[15].m13,
            ]),
            m22: simba::simd::f32x16::from([
                data[0].m22,
                data[1].m22,
                data[2].m22,
                data[3].m22,
                data[4].m22,
                data[5].m22,
                data[6].m22,
                data[7].m22,
                data[8].m22,
                data[9].m22,
                data[10].m22,
                data[11].m22,
                data[12].m22,
                data[13].m22,
                data[14].m22,
                data[15].m22,
            ]),
            m23: simba::simd::f32x16::from([
                data[0].m23,
                data[1].m23,
                data[2].m23,
                data[3].m23,
                data[4].m23,
                data[5].m23,
                data[6].m23,
                data[7].m23,
                data[8].m23,
                data[9].m23,
                data[10].m23,
                data[11].m23,
                data[12].m23,
                data[13].m23,
                data[14].m23,
                data[15].m23,
            ]),
            m33: simba::simd::f32x16::from([
                data[0].m33,
                data[1].m33,
                data[2].m33,
                data[3].m33,
                data[4].m33,
                data[5].m33,
                data[6].m33,
                data[7].m33,
                data[8].m33,
                data[9].m33,
                data[10].m33,
                data[11].m33,
                data[12].m33,
                data[13].m33,
                data[14].m33,
                data[15].m33,
            ]),
        }
    }
}
