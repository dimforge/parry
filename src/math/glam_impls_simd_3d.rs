use core::ops::{
    Add, AddAssign, Div, DivAssign, Index, IndexMut, Mul, MulAssign, Neg, Rem, RemAssign, Sub,
    SubAssign,
};

use crate::math::{Iso3, PointOps, Real, SimdBool, SimdReal, SIMD_WIDTH};
use crate::utils::SdpMatrix3;
use glam::{Quat, Vec3};
use simba::simd::{SimdComplexField, SimdPartialOrd, SimdValue};

#[derive(Copy, Clone, Debug)]
pub struct SimdMat3 {
    pub x_axis: SimdVec3,
    pub y_axis: SimdVec3,
    pub z_axis: SimdVec3,
}

impl SimdMat3 {
    #[inline]
    pub fn column(&self, i: usize) -> &SimdVec3 {
        [&self.x_axis, &self.y_axis, &self.z_axis][i]
    }

    #[inline]
    pub fn transpose(&self) -> Self {
        Self {
            x_axis: SimdVec3::new(self.x_axis.x, self.y_axis.x, self.z_axis.x),
            y_axis: SimdVec3::new(self.x_axis.y, self.y_axis.y, self.z_axis.y),
            z_axis: SimdVec3::new(self.x_axis.z, self.y_axis.z, self.z_axis.z),
        }
    }
}

impl Index<(usize, usize)> for SimdMat3 {
    type Output = SimdReal;

    #[inline]
    fn index(&self, (irow, icol): (usize, usize)) -> &Self::Output {
        &[&self.x_axis, &self.y_axis, &self.z_axis][icol][irow]
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SimdVec3 {
    pub x: SimdReal,
    pub y: SimdReal,
    pub z: SimdReal,
}

impl Default for SimdVec3 {
    #[inline]
    fn default() -> Self {
        Self::zeros()
    }
}

impl Index<usize> for SimdVec3 {
    type Output = SimdReal;
    #[inline]
    fn index(&self, index: usize) -> &Self::Output {
        [&self.x, &self.y, &self.z][index]
    }
}

impl IndexMut<usize> for SimdVec3 {
    #[inline]
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        [&mut self.x, &mut self.y, &mut self.z][index]
    }
}

impl From<[Vec3; SIMD_WIDTH]> for SimdVec3 {
    #[inline]
    fn from(value: [Vec3; SIMD_WIDTH]) -> Self {
        Self {
            x: array![|ii| value[ii].x; SIMD_WIDTH].into(),
            y: array![|ii| value[ii].y; SIMD_WIDTH].into(),
            z: array![|ii| value[ii].z; SIMD_WIDTH].into(),
        }
    }
}

impl From<[SimdReal; 3]> for SimdVec3 {
    #[inline]
    fn from(value: [SimdReal; 3]) -> Self {
        Self {
            x: value[0],
            y: value[1],
            z: value[2],
        }
    }
}

impl From<SimdVec3> for [SimdReal; 3] {
    #[inline]
    fn from(value: SimdVec3) -> Self {
        [value.x, value.y, value.z]
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SimdQuat {
    pub x: SimdReal,
    pub y: SimdReal,
    pub z: SimdReal,
    pub w: SimdReal,
}

impl Default for SimdQuat {
    #[inline]
    fn default() -> Self {
        Self::identity()
    }
}

impl From<[Quat; 4]> for SimdQuat {
    fn from(value: [Quat; 4]) -> Self {
        Self {
            x: value.map(|q| q.as_ref()[0]).into(),
            y: value.map(|q| q.as_ref()[1]).into(),
            z: value.map(|q| q.as_ref()[2]).into(),
            w: value.map(|q| q.as_ref()[3]).into(),
        }
    }
}

impl From<SimdQuat> for SimdMat3 {
    #[inline]
    #[must_use]
    fn from(rotation: SimdQuat) -> Self {
        let x2 = rotation.x + rotation.x;
        let y2 = rotation.y + rotation.y;
        let z2 = rotation.z + rotation.z;
        let xx = rotation.x * x2;
        let xy = rotation.x * y2;
        let xz = rotation.x * z2;
        let yy = rotation.y * y2;
        let yz = rotation.y * z2;
        let zz = rotation.z * z2;
        let wx = rotation.w * x2;
        let wy = rotation.w * y2;
        let wz = rotation.w * z2;

        let one = SimdReal::splat(1.0);

        Self {
            x_axis: SimdVec3 {
                x: one - (yy + zz),
                y: xy + wz,
                z: xz - wy,
            },
            y_axis: SimdVec3 {
                x: xy - wz,
                y: one - (xx + zz),
                z: yz + wx,
            },
            z_axis: SimdVec3 {
                x: xz + wy,
                y: yz - wx,
                z: one - (xx + yy),
            },
        }
    }
}

impl Mul<SimdVec3> for SimdMat3 {
    type Output = SimdVec3;
    #[inline]
    fn mul(self, rhs: SimdVec3) -> Self::Output {
        self.x_axis * rhs.x + self.y_axis * rhs.y + self.z_axis * rhs.z
    }
}

impl Mul<SimdMat3> for SimdMat3 {
    type Output = Self;
    #[inline]
    fn mul(self, rhs: SimdMat3) -> Self::Output {
        Self {
            x_axis: self.x_axis * rhs.x_axis.x
                + self.y_axis * rhs.x_axis.y
                + self.z_axis * rhs.x_axis.z,
            y_axis: self.x_axis * rhs.y_axis.x
                + self.y_axis * rhs.y_axis.y
                + self.z_axis * rhs.y_axis.z,
            z_axis: self.x_axis * rhs.z_axis.x
                + self.y_axis * rhs.z_axis.y
                + self.z_axis * rhs.z_axis.z,
        }
    }
}

impl Mul<SimdReal> for SimdMat3 {
    type Output = Self;
    #[inline]
    fn mul(self, rhs: SimdReal) -> Self::Output {
        Self {
            x_axis: self.x_axis * rhs,
            y_axis: self.y_axis * rhs,
            z_axis: self.z_axis * rhs,
        }
    }
}

impl MulAssign<SimdReal> for SimdMat3 {
    #[inline]
    fn mul_assign(&mut self, rhs: SimdReal) {
        self.x_axis *= rhs;
        self.y_axis *= rhs;
        self.z_axis *= rhs;
    }
}

impl Add<SimdMat3> for SimdMat3 {
    type Output = Self;
    #[inline]
    fn add(self, rhs: SimdMat3) -> Self::Output {
        Self {
            x_axis: self.x_axis + rhs.x_axis,
            y_axis: self.y_axis + rhs.y_axis,
            z_axis: self.z_axis + rhs.z_axis,
        }
    }
}

impl Sub<SimdMat3> for SimdMat3 {
    type Output = Self;
    #[inline]
    fn sub(self, rhs: SimdMat3) -> Self::Output {
        Self {
            x_axis: self.x_axis - rhs.x_axis,
            y_axis: self.y_axis - rhs.y_axis,
            z_axis: self.z_axis - rhs.z_axis,
        }
    }
}

impl SimdMat3 {
    #[inline]
    #[rustfmt::skip]
    pub fn new(
        m11: SimdReal, m12: SimdReal, m13: SimdReal,
        m21: SimdReal, m22: SimdReal, m23: SimdReal,
        m31: SimdReal, m32: SimdReal, m33: SimdReal,
    ) -> Self {
        Self {
            x_axis: SimdVec3::new(m11, m21, m31),
            y_axis: SimdVec3::new(m12, m22, m32),
            z_axis: SimdVec3::new(m13, m23, m33),
        }
    }

    #[inline]
    #[rustfmt::skip]
    pub fn from_diagonal_element(elt: SimdReal) -> Self {
        let zero = SimdReal::splat(0.0);
        Self::new(
            elt, zero, zero,
            zero, elt, zero,
            zero, zero, elt
        )
    }

    #[inline]
    pub fn abs(self) -> Self {
        Self {
            x_axis: self.x_axis.abs(),
            y_axis: self.y_axis.abs(),
            z_axis: self.z_axis.abs(),
        }
    }

    #[inline]
    pub fn map(self, mut f: impl FnMut(SimdReal) -> SimdReal) -> Self {
        Self {
            x_axis: SimdVec3 {
                x: f(self.x_axis.x),
                y: f(self.x_axis.y),
                z: f(self.x_axis.z),
            },
            y_axis: SimdVec3 {
                x: f(self.y_axis.x),
                y: f(self.y_axis.y),
                z: f(self.y_axis.z),
            },
            z_axis: SimdVec3 {
                x: f(self.z_axis.x),
                y: f(self.z_axis.y),
                z: f(self.z_axis.z),
            },
        }
    }

    #[inline]
    #[doc(hidden)]
    pub fn into_inner(self) -> Self {
        self
    }
}

impl SimdVec3 {
    #[inline]
    pub fn new(x: SimdReal, y: SimdReal, z: SimdReal) -> Self {
        Self { x, y, z }
    }

    #[inline]
    pub fn splat(vec: Vec3) -> Self {
        Self {
            x: SimdReal::splat(vec.x),
            y: SimdReal::splat(vec.y),
            z: SimdReal::splat(vec.z),
        }
    }

    #[inline]
    pub const fn repeat(val: SimdReal) -> Self {
        Self {
            x: val,
            y: val,
            z: val,
        }
    }

    #[inline]
    pub const fn splat_simd_real(real: SimdReal) -> Self {
        Self::repeat(real)
    }

    #[inline]
    pub fn zeros() -> Self {
        Self::repeat(SimdReal::splat(0.0))
    }

    #[inline]
    pub fn from_vecs(vecs: [Vec3; 4]) -> Self {
        Self {
            x: SimdReal::from([vecs[0].x, vecs[1].x, vecs[2].x, vecs[3].x]),
            y: SimdReal::from([vecs[0].y, vecs[1].y, vecs[2].y, vecs[3].y]),
            z: SimdReal::from([vecs[0].z, vecs[1].z, vecs[2].z, vecs[3].z]),
        }
    }

    #[inline]
    pub fn extract(&self, lane: usize) -> Vec3 {
        Vec3 {
            x: self.x.extract(lane),
            y: self.y.extract(lane),
            z: self.z.extract(lane),
        }
    }

    #[inline]
    pub fn map<F: Fn(SimdReal) -> f32>(self, f: F) -> Vec3 {
        Vec3 {
            x: f(self.x),
            y: f(self.y),
            z: f(self.z),
        }
    }

    #[inline]
    pub fn replace(&mut self, lane: usize, new: Vec3) {
        self.x.replace(lane, new.x);
        self.y.replace(lane, new.y);
        self.z.replace(lane, new.z);
    }

    #[inline]
    pub fn select(self, cond: SimdBool, other: Self) -> Self {
        Self {
            x: self.x.select(cond, other.x),
            y: self.y.select(cond, other.y),
            z: self.z.select(cond, other.z),
        }
    }

    #[inline]
    pub fn length(self) -> SimdReal {
        self.length_squared().simd_sqrt()
    }

    #[inline]
    pub fn length_squared(self) -> SimdReal {
        self.dot(self)
    }

    #[inline]
    pub fn dot(self, other: Self) -> SimdReal {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    #[inline]
    pub fn cross(self, other: Self) -> SimdVec3 {
        Self {
            x: self.y * other.z - self.z * other.y,
            y: self.z * other.x - self.x * other.z,
            z: self.x * other.y - self.y * other.x,
        }
    }

    #[inline]
    pub fn min(self, other: Self) -> Self {
        Self {
            x: self.x.simd_min(other.x),
            y: self.y.simd_min(other.y),
            z: self.y.simd_min(other.z),
        }
    }

    #[inline]
    pub fn max(self, other: Self) -> Self {
        Self {
            x: self.x.simd_max(other.x),
            y: self.y.simd_max(other.y),
            z: self.y.simd_max(other.z),
        }
    }

    #[inline]
    pub fn abs(self) -> Self {
        Self {
            x: self.x.simd_abs(),
            y: self.y.simd_abs(),
            z: self.z.simd_abs(),
        }
    }

    #[inline]
    pub fn outer_product(self, rhs: Self) -> SimdMat3 {
        SimdMat3 {
            x_axis: self * rhs.x,
            y_axis: self * rhs.y,
            z_axis: self * rhs.z,
        }
    }

    /*
     * The functions below are just useful for the glam/nalgebra api compatibility.
     */
    #[inline]
    #[doc(hidden)]
    pub fn into_owned(self) -> Self {
        self
    }

    #[inline]
    #[doc(hidden)]
    pub fn coords(self) -> Self {
        self
    }

    #[inline]
    #[doc(hidden)]
    pub fn norm(&self) -> SimdReal {
        self.length()
    }

    #[inline]
    #[doc(hidden)]
    pub fn normalize_mut(&mut self) -> SimdReal {
        let norm = self.norm();
        *self /= norm;
        norm
    }

    #[inline]
    #[doc(hidden)]
    pub fn inf(&self, rhs: &Self) -> Self {
        self.min(*rhs)
    }

    #[inline]
    #[doc(hidden)]
    pub fn sup(&self, rhs: &Self) -> Self {
        self.max(*rhs)
    }

    #[inline]
    #[doc(hidden)]
    pub fn component_mul(self, rhs: &Self) -> Self {
        self * *rhs
    }
}

impl PointOps for SimdVec3 {
    type Coords = Self;
    #[inline(always)]
    fn into_vector(self) -> Self::Coords {
        self
    }
    #[inline(always)]
    fn as_vector_mut(&mut self) -> &mut Self::Coords {
        self
    }
    #[inline(always)]
    fn origin() -> Self {
        Self::zeros()
    }
}

impl Add<SimdVec3> for SimdVec3 {
    type Output = Self;
    #[inline]
    fn add(self, rhs: Self) -> Self {
        Self {
            x: self.x.add(rhs.x),
            y: self.y.add(rhs.y),
            z: self.z.add(rhs.z),
        }
    }
}

impl AddAssign<SimdVec3> for SimdVec3 {
    #[inline]
    fn add_assign(&mut self, rhs: Self) {
        self.x.add_assign(rhs.x);
        self.y.add_assign(rhs.y);
        self.z.add_assign(rhs.z);
    }
}

impl Add<SimdReal> for SimdVec3 {
    type Output = Self;
    #[inline]
    fn add(self, rhs: SimdReal) -> Self {
        Self {
            x: self.x.add(rhs),
            y: self.y.add(rhs),
            z: self.z.add(rhs),
        }
    }
}

impl AddAssign<SimdReal> for SimdVec3 {
    #[inline]
    fn add_assign(&mut self, rhs: SimdReal) {
        self.x.add_assign(rhs);
        self.y.add_assign(rhs);
        self.z.add_assign(rhs);
    }
}

impl Add<SimdVec3> for SimdReal {
    type Output = SimdVec3;
    #[inline]
    fn add(self, rhs: SimdVec3) -> SimdVec3 {
        SimdVec3 {
            x: self.add(rhs.x),
            y: self.add(rhs.y),
            z: self.add(rhs.z),
        }
    }
}

impl Sub<SimdVec3> for SimdVec3 {
    type Output = Self;
    #[inline]
    fn sub(self, rhs: Self) -> Self {
        Self {
            x: self.x.sub(rhs.x),
            y: self.y.sub(rhs.y),
            z: self.z.sub(rhs.z),
        }
    }
}

impl SubAssign<SimdVec3> for SimdVec3 {
    #[inline]
    fn sub_assign(&mut self, rhs: SimdVec3) {
        self.x.sub_assign(rhs.x);
        self.y.sub_assign(rhs.y);
        self.z.sub_assign(rhs.z);
    }
}

impl Sub<SimdReal> for SimdVec3 {
    type Output = Self;
    #[inline]
    fn sub(self, rhs: SimdReal) -> Self {
        Self {
            x: self.x.sub(rhs),
            y: self.y.sub(rhs),
            z: self.z.sub(rhs),
        }
    }
}

impl SubAssign<SimdReal> for SimdVec3 {
    #[inline]
    fn sub_assign(&mut self, rhs: SimdReal) {
        self.x.sub_assign(rhs);
        self.y.sub_assign(rhs);
        self.z.sub_assign(rhs);
    }
}

impl Sub<SimdVec3> for SimdReal {
    type Output = SimdVec3;
    #[inline]
    fn sub(self, rhs: SimdVec3) -> SimdVec3 {
        SimdVec3 {
            x: self.sub(rhs.x),
            y: self.sub(rhs.y),
            z: self.sub(rhs.z),
        }
    }
}

impl Mul<SimdVec3> for SimdVec3 {
    type Output = Self;
    #[inline]
    fn mul(self, rhs: Self) -> Self {
        Self {
            x: self.x.mul(rhs.x),
            y: self.y.mul(rhs.y),
            z: self.z.mul(rhs.z),
        }
    }
}

impl MulAssign<SimdVec3> for SimdVec3 {
    #[inline]
    fn mul_assign(&mut self, rhs: Self) {
        self.x.mul_assign(rhs.x);
        self.y.mul_assign(rhs.y);
        self.z.mul_assign(rhs.z);
    }
}

impl<'a> Mul<SimdReal> for &'a SimdVec3 {
    type Output = SimdVec3;
    #[inline]
    fn mul(self, rhs: SimdReal) -> SimdVec3 {
        SimdVec3 {
            x: self.x.mul(rhs),
            y: self.y.mul(rhs),
            z: self.z.mul(rhs),
        }
    }
}

impl Mul<SimdReal> for SimdVec3 {
    type Output = Self;
    #[inline]
    fn mul(self, rhs: SimdReal) -> Self {
        Self {
            x: self.x.mul(rhs),
            y: self.y.mul(rhs),
            z: self.z.mul(rhs),
        }
    }
}

impl MulAssign<SimdReal> for SimdVec3 {
    #[inline]
    fn mul_assign(&mut self, rhs: SimdReal) {
        self.x.mul_assign(rhs);
        self.y.mul_assign(rhs);
        self.z.mul_assign(rhs);
    }
}

impl Mul<SimdVec3> for SimdReal {
    type Output = SimdVec3;
    #[inline]
    fn mul(self, rhs: SimdVec3) -> SimdVec3 {
        SimdVec3 {
            x: self.mul(rhs.x),
            y: self.mul(rhs.y),
            z: self.mul(rhs.z),
        }
    }
}

impl Div<SimdVec3> for SimdVec3 {
    type Output = Self;
    #[inline]
    fn div(self, rhs: Self) -> Self {
        Self {
            x: self.x.div(rhs.x),
            y: self.y.div(rhs.y),
            z: self.z.div(rhs.z),
        }
    }
}

impl DivAssign<SimdVec3> for SimdVec3 {
    #[inline]
    fn div_assign(&mut self, rhs: Self) {
        self.x.div_assign(rhs.x);
        self.y.div_assign(rhs.y);
        self.z.div_assign(rhs.z);
    }
}

impl Div<SimdReal> for SimdVec3 {
    type Output = Self;
    #[inline]
    fn div(self, rhs: SimdReal) -> Self {
        Self {
            x: self.x.div(rhs),
            y: self.y.div(rhs),
            z: self.z.div(rhs),
        }
    }
}

impl DivAssign<SimdReal> for SimdVec3 {
    #[inline]
    fn div_assign(&mut self, rhs: SimdReal) {
        self.x.div_assign(rhs);
        self.y.div_assign(rhs);
        self.z.div_assign(rhs);
    }
}

impl Div<SimdVec3> for SimdReal {
    type Output = SimdVec3;
    #[inline]
    fn div(self, rhs: SimdVec3) -> SimdVec3 {
        SimdVec3 {
            x: self.div(rhs.x),
            y: self.div(rhs.y),
            z: self.div(rhs.z),
        }
    }
}

impl Rem<SimdVec3> for SimdVec3 {
    type Output = Self;
    #[inline]
    fn rem(self, rhs: Self) -> Self {
        Self {
            x: self.x.rem(rhs.x),
            y: self.y.rem(rhs.y),
            z: self.z.rem(rhs.z),
        }
    }
}

impl RemAssign<SimdVec3> for SimdVec3 {
    #[inline]
    fn rem_assign(&mut self, rhs: Self) {
        self.x.rem_assign(rhs.x);
        self.y.rem_assign(rhs.y);
        self.z.rem_assign(rhs.z);
    }
}

impl Rem<SimdReal> for SimdVec3 {
    type Output = Self;
    #[inline]
    fn rem(self, rhs: SimdReal) -> Self {
        Self {
            x: self.x.rem(rhs),
            y: self.y.rem(rhs),
            z: self.z.rem(rhs),
        }
    }
}

impl RemAssign<SimdReal> for SimdVec3 {
    #[inline]
    fn rem_assign(&mut self, rhs: SimdReal) {
        self.x.rem_assign(rhs);
        self.y.rem_assign(rhs);
        self.z.rem_assign(rhs);
    }
}

impl Rem<SimdVec3> for SimdReal {
    type Output = SimdVec3;
    #[inline]
    fn rem(self, rhs: SimdVec3) -> SimdVec3 {
        SimdVec3 {
            x: self.rem(rhs.x),
            y: self.rem(rhs.y),
            z: self.rem(rhs.z),
        }
    }
}

impl Neg for SimdVec3 {
    type Output = Self;
    #[inline]
    fn neg(self) -> Self::Output {
        Self {
            x: -self.x,
            y: -self.y,
            z: -self.z,
        }
    }
}

impl SimdQuat {
    #[inline]
    pub fn splat(quat: Quat) -> Self {
        Self {
            x: SimdReal::splat(quat.x),
            y: SimdReal::splat(quat.y),
            z: SimdReal::splat(quat.z),
            w: SimdReal::splat(quat.w),
        }
    }

    #[inline]
    pub fn identity() -> Self {
        Self::splat(Quat::IDENTITY)
    }

    #[inline]
    pub fn extract(&self, lane: usize) -> Quat {
        Quat::from_xyzw(
            self.x.extract(lane),
            self.y.extract(lane),
            self.z.extract(lane),
            self.w.extract(lane),
        )
    }

    #[inline]
    pub fn from_scaled_axis(axis: Vec3) -> Self {
        Self::splat(Quat::from_scaled_axis(axis))
    }

    #[inline]
    pub fn from_axis_angle(axis: Vec3, angle: Real) -> Self {
        Self::splat(Quat::from_axis_angle(axis, angle))
    }

    #[inline]
    pub fn imag(&self) -> SimdVec3 {
        SimdVec3::new(self.x, self.y, self.z)
    }

    #[inline]
    pub fn dot(&self, rhs: Self) -> SimdReal {
        self.x * rhs.x + self.y * rhs.y + self.z * rhs.z + self.w * rhs.w
    }

    #[inline]
    pub fn normalize(self) -> Self {
        let length = self.length();
        Self {
            x: self.x / length,
            y: self.y / length,
            z: self.z / length,
            w: self.w / length,
        }
    }

    #[inline]
    pub fn length(&self) -> SimdReal {
        (self.x * self.x + self.y * self.y + self.y * self.y + self.z * self.z).simd_sqrt()
    }

    #[inline]
    pub fn inverse(self) -> Self {
        Self {
            x: -self.x,
            y: -self.y,
            z: -self.z,
            w: self.w,
        }
    }

    #[inline]
    pub fn to_rotation_matrix(self) -> SimdMat3 {
        self.into()
    }

    #[inline]
    pub fn xyz(&self) -> SimdVec3 {
        SimdVec3::new(self.x, self.y, self.z)
    }

    #[doc(hidden)]
    #[inline]
    pub fn into_inner(self) -> Self {
        self
    }

    #[doc(hidden)]
    #[inline]
    pub fn as_mut_unchecked(&mut self) -> &mut Self {
        self
    }
}

impl Add<SimdQuat> for SimdQuat {
    type Output = Self;
    #[inline]
    fn add(self, rhs: Self) -> Self {
        Self {
            x: self.x + rhs.x,
            y: self.y + rhs.y,
            z: self.z + rhs.z,
            w: self.w + rhs.w,
        }
    }
}

impl Sub<SimdQuat> for SimdQuat {
    type Output = Self;
    #[inline]
    fn sub(self, rhs: Self) -> Self {
        Self {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
            z: self.z - rhs.z,
            w: self.w - rhs.w,
        }
    }
}

impl Mul<SimdReal> for SimdQuat {
    type Output = Self;
    #[inline]
    fn mul(self, rhs: SimdReal) -> Self {
        Self {
            x: self.x * rhs,
            y: self.y * rhs,
            z: self.z * rhs,
            w: self.w * rhs,
        }
    }
}

impl Div<SimdReal> for SimdQuat {
    type Output = Self;
    #[inline]
    fn div(self, rhs: SimdReal) -> Self {
        Self {
            x: self.x / rhs,
            y: self.y / rhs,
            z: self.z / rhs,
            w: self.w / rhs,
        }
    }
}

impl Mul<SimdQuat> for SimdQuat {
    type Output = SimdQuat;
    #[inline]
    fn mul(self, rhs: SimdQuat) -> Self::Output {
        Self {
            x: self.w * rhs.x + self.x * rhs.w + self.y * rhs.z - self.z * rhs.y,
            y: self.w * rhs.y - self.x * rhs.z + self.y * rhs.w + self.z * rhs.x,
            z: self.w * rhs.z + self.x * rhs.y - self.y * rhs.x + self.z * rhs.w,
            w: self.w * rhs.w - self.x * rhs.x - self.y * rhs.y - self.z * rhs.z,
        }
    }
}

impl Mul<SimdVec3> for SimdQuat {
    type Output = SimdVec3;
    #[inline]
    fn mul(self, v: SimdVec3) -> Self::Output {
        let u = SimdVec3 {
            x: self.x,
            y: self.y,
            z: self.z,
        };
        let two = SimdReal::splat(2.0);
        two * u.dot(v) * u + (self.w * self.w - u.dot(u)) * v + two * self.w * u.cross(v)
    }
}

impl Neg for SimdQuat {
    type Output = Self;
    #[inline]
    fn neg(self) -> Self {
        self * SimdReal::splat(-1.0)
    }
}

#[derive(Clone, Copy, Default, Debug, PartialEq)]
pub struct SimdIso3 {
    pub translation: SimdVec3,
    pub rotation: SimdQuat,
}

impl SimdIso3 {
    // /// An identity isometry.
    // pub const IDENTITY: Self = Self {
    //     translation: SimdVec3::zeros(),
    //     rotation: SimdQuat::IDENTITY,
    // };

    #[inline]
    pub fn new(translation: Vec3, rotation: Vec3) -> Self {
        Self {
            translation: SimdVec3::splat(translation),
            rotation: SimdQuat::from_scaled_axis(rotation),
        }
    }

    #[inline]
    pub fn identity() -> Self {
        Self {
            translation: SimdVec3::zeros(),
            rotation: SimdQuat::identity(),
        }
    }

    #[inline]
    pub fn splat(iso: Iso3) -> Self {
        Self {
            translation: SimdVec3::splat(iso.translation),
            rotation: SimdQuat::splat(iso.rotation),
        }
    }

    #[inline]
    pub fn extract(&self, lane: usize) -> Iso3 {
        Iso3 {
            translation: self.translation.extract(lane),
            rotation: self.rotation.extract(lane),
        }
    }

    #[inline]
    pub fn from_translation(translation: Vec3) -> Self {
        Self {
            translation: SimdVec3::splat(translation),
            rotation: SimdQuat::identity(),
        }
    }

    #[inline]
    pub fn from_rotation(rotation: Quat) -> Self {
        Self {
            translation: SimdVec3::zeros(),
            rotation: SimdQuat::splat(rotation),
        }
    }

    #[inline]
    pub fn inverse(self) -> Self {
        let inv_rot = self.rotation.inverse();
        Self {
            translation: inv_rot * -self.translation,
            rotation: inv_rot,
        }
    }

    #[inline]
    pub fn transform_point(self, point: SimdVec3) -> SimdVec3 {
        self.translation + self.rotation * point
    }

    #[inline]
    pub fn inverse_transform_point(self, point: SimdVec3) -> SimdVec3 {
        self.rotation.inverse() * (point - self.translation)
    }
}

impl From<[Iso3; 4]> for SimdIso3 {
    fn from(value: [Iso3; 4]) -> Self {
        Self {
            rotation: value.map(|iso| iso.rotation).into(),
            translation: value.map(|tra| tra.translation).into(),
        }
    }
}

impl Mul<SimdIso3> for SimdIso3 {
    type Output = SimdIso3;
    #[inline]
    fn mul(self, rhs: SimdIso3) -> Self::Output {
        Self {
            rotation: self.rotation * rhs.rotation,
            translation: self.translation + self.rotation * rhs.translation,
        }
    }
}

impl Mul<SimdVec3> for SdpMatrix3<SimdReal> {
    type Output = SimdVec3;

    #[inline]
    fn mul(self, rhs: SimdVec3) -> Self::Output {
        let x = self.m11 * rhs.x + self.m12 * rhs.y + self.m13 * rhs.z;
        let y = self.m12 * rhs.x + self.m22 * rhs.y + self.m23 * rhs.z;
        let z = self.m13 * rhs.x + self.m23 * rhs.y + self.m33 * rhs.z;
        SimdVec3::new(x, y, z)
    }
}
