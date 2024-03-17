use core::ops::{
    Add, AddAssign, Div, DivAssign, Index, IndexMut, Mul, MulAssign, Neg, Rem, RemAssign, Sub,
    SubAssign,
};

use crate::math::{IntoInner, Iso2, PointOps, Real, SimdBool, SimdReal, SIMD_WIDTH};
use crate::utils::SdpMatrix2;
use glam::{Mat2, Vec2};
use num::{One, Zero};
use simba::simd::{SimdComplexField, SimdPartialOrd, SimdValue};

#[derive(Copy, Clone, Debug, PartialEq)]
pub struct SimdMat2 {
    pub x_axis: SimdVec2,
    pub y_axis: SimdVec2,
}

impl Default for SimdMat2 {
    fn default() -> Self {
        Self::identity()
    }
}

impl SimdMat2 {
    #[inline]
    pub fn identity() -> Self {
        let zero = SimdReal::zero();
        let one = SimdReal::one();

        Self {
            x_axis: SimdVec2::new(one, zero),
            y_axis: SimdVec2::new(zero, one),
        }
    }

    #[inline]
    pub fn from_angle(ang: f32) -> Self {
        let simd_ang = SimdReal::splat(ang);
        let (s, c) = simd_ang.simd_sin_cos();

        SimdMat2 {
            x_axis: SimdVec2::new(c, s),
            y_axis: SimdVec2::new(-s, c),
        }
    }

    #[inline]
    pub fn splat(mat2: Mat2) -> Self {
        SimdMat2 {
            x_axis: SimdVec2::splat(mat2.col(0)),
            y_axis: SimdVec2::splat(mat2.col(1)),
        }
    }

    #[inline]
    pub fn column(&self, i: usize) -> &SimdVec2 {
        [&self.x_axis, &self.y_axis][i]
    }

    #[inline]
    pub fn transpose(&self) -> Self {
        Self {
            x_axis: SimdVec2::new(self.x_axis.x, self.y_axis.x),
            y_axis: SimdVec2::new(self.x_axis.y, self.y_axis.y),
        }
    }

    #[inline]
    pub fn extract(&self, i: usize) -> Mat2 {
        Mat2::from_cols(self.x_axis.extract(i), self.y_axis.extract(i))
    }

    #[doc(hidden)]
    #[inline]
    pub fn to_rotation_matrix(self) -> Self {
        self
    }

    #[inline]
    pub fn imag(&self) -> SimdReal {
        self.x_axis.y
    }
}

impl Index<(usize, usize)> for SimdMat2 {
    type Output = SimdReal;

    #[inline]
    fn index(&self, (irow, icol): (usize, usize)) -> &Self::Output {
        &[&self.x_axis, &self.y_axis][icol][irow]
    }
}

impl crate::math::Mat2Ops for SimdMat2 {
    type Element = SimdReal;
    #[inline]
    fn angle(&self) -> SimdReal {
        use na::SimdRealField;
        self.x_axis.y.clone().simd_atan2(self.x_axis.x.clone())
    }
}

impl From<[Mat2; 4]> for SimdMat2 {
    #[inline]
    fn from(value: [Mat2; 4]) -> Self {
        Self {
            x_axis: value.map(|m| m.x_axis).into(),
            y_axis: value.map(|m| m.y_axis).into(),
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Default)]
pub struct SimdVec2 {
    pub x: SimdReal,
    pub y: SimdReal,
}

impl Index<usize> for SimdVec2 {
    type Output = SimdReal;
    #[inline]
    fn index(&self, index: usize) -> &Self::Output {
        [&self.x, &self.y][index]
    }
}

impl IndexMut<usize> for SimdVec2 {
    #[inline]
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        [&mut self.x, &mut self.y][index]
    }
}

impl From<[Vec2; SIMD_WIDTH]> for SimdVec2 {
    #[inline]
    fn from(value: [Vec2; SIMD_WIDTH]) -> Self {
        Self {
            x: array![|ii| value[ii].x; SIMD_WIDTH].into(),
            y: array![|ii| value[ii].y; SIMD_WIDTH].into(),
        }
    }
}

impl From<[SimdReal; 2]> for SimdVec2 {
    #[inline]
    fn from(value: [SimdReal; 2]) -> Self {
        Self {
            x: value[0],
            y: value[1],
        }
    }
}

impl From<SimdVec2> for [SimdReal; 2] {
    #[inline]
    fn from(value: SimdVec2) -> Self {
        [value.x, value.y]
    }
}

impl Mul<SimdVec2> for SimdMat2 {
    type Output = SimdVec2;
    #[inline]
    fn mul(self, rhs: SimdVec2) -> Self::Output {
        self.x_axis * rhs.x + self.y_axis * rhs.y
    }
}

impl Mul<SimdMat2> for SimdMat2 {
    type Output = Self;
    #[inline]
    fn mul(self, rhs: SimdMat2) -> Self::Output {
        Self {
            x_axis: self.x_axis * rhs.x_axis.x + self.y_axis * rhs.x_axis.y,
            y_axis: self.x_axis * rhs.y_axis.x + self.y_axis * rhs.y_axis.y,
        }
    }
}

impl Mul<SimdReal> for SimdMat2 {
    type Output = Self;
    #[inline]
    fn mul(self, rhs: SimdReal) -> Self::Output {
        Self {
            x_axis: self.x_axis * rhs,
            y_axis: self.y_axis * rhs,
        }
    }
}

impl MulAssign<SimdReal> for SimdMat2 {
    #[inline]
    fn mul_assign(&mut self, rhs: SimdReal) {
        self.x_axis *= rhs;
        self.y_axis *= rhs;
    }
}

impl Add<SimdMat2> for SimdMat2 {
    type Output = Self;
    #[inline]
    fn add(self, rhs: SimdMat2) -> Self::Output {
        Self {
            x_axis: self.x_axis + rhs.x_axis,
            y_axis: self.y_axis + rhs.y_axis,
        }
    }
}

impl Sub<SimdMat2> for SimdMat2 {
    type Output = Self;
    #[inline]
    fn sub(self, rhs: SimdMat2) -> Self::Output {
        Self {
            x_axis: self.x_axis - rhs.x_axis,
            y_axis: self.y_axis - rhs.y_axis,
        }
    }
}

impl SimdMat2 {
    #[inline]
    #[rustfmt::skip]
    pub fn new(
        m11: SimdReal, m12: SimdReal,
        m21: SimdReal, m22: SimdReal,
    ) -> Self {
        Self {
            x_axis: SimdVec2::new(m11, m21),
            y_axis: SimdVec2::new(m12, m22),
        }
    }

    #[inline]
    #[rustfmt::skip]
    pub fn from_diagonal_element(elt: SimdReal) -> Self {
        let zero = SimdReal::splat(0.0);
        Self::new(
            elt, zero,
            zero, elt,
        )
    }

    #[inline]
    pub fn abs(self) -> Self {
        Self {
            x_axis: self.x_axis.abs(),
            y_axis: self.y_axis.abs(),
        }
    }

    #[inline]
    pub fn map(self, mut f: impl FnMut(SimdReal) -> SimdReal) -> Self {
        Self {
            x_axis: SimdVec2 {
                x: f(self.x_axis.x),
                y: f(self.x_axis.y),
            },
            y_axis: SimdVec2 {
                x: f(self.y_axis.x),
                y: f(self.y_axis.y),
            },
        }
    }

    #[inline]
    #[doc(hidden)]
    pub fn into_inner(self) -> Self {
        self
    }
}

impl SimdVec2 {
    #[inline]
    pub fn new(x: SimdReal, y: SimdReal) -> Self {
        Self { x, y }
    }

    #[inline]
    pub fn splat(vec: Vec2) -> Self {
        Self {
            x: SimdReal::splat(vec.x),
            y: SimdReal::splat(vec.y),
        }
    }

    #[inline]
    pub const fn repeat(val: SimdReal) -> Self {
        Self { x: val, y: val }
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
    pub fn from_vecs(vecs: [Vec2; 4]) -> Self {
        Self {
            x: SimdReal::from([vecs[0].x, vecs[1].x, vecs[2].x, vecs[2].x]),
            y: SimdReal::from([vecs[0].y, vecs[1].y, vecs[2].y, vecs[2].y]),
        }
    }

    #[inline]
    pub fn extract(&self, lane: usize) -> Vec2 {
        Vec2 {
            x: self.x.extract(lane),
            y: self.y.extract(lane),
        }
    }

    #[inline]
    pub fn map<F: Fn(SimdReal) -> f32>(self, f: F) -> Vec2 {
        Vec2 {
            x: f(self.x),
            y: f(self.y),
        }
    }

    #[inline]
    pub fn replace(&mut self, lane: usize, new: Vec2) {
        self.x.replace(lane, new.x);
        self.y.replace(lane, new.y);
    }

    #[inline]
    pub fn select(self, cond: SimdBool, other: Self) -> Self {
        Self {
            x: self.x.select(cond, other.x),
            y: self.y.select(cond, other.y),
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
        self.x * other.x + self.y * other.y
    }

    #[inline]
    pub fn cross(self, other: Self) -> SimdReal {
        self.x * other.y - self.y * other.x
    }

    #[inline]
    pub fn min(self, other: Self) -> Self {
        Self {
            x: self.x.simd_min(other.x),
            y: self.y.simd_min(other.y),
        }
    }

    #[inline]
    pub fn max(self, other: Self) -> Self {
        Self {
            x: self.x.simd_max(other.x),
            y: self.y.simd_max(other.y),
        }
    }

    #[inline]
    pub fn abs(self) -> Self {
        Self {
            x: self.x.simd_abs(),
            y: self.y.simd_abs(),
        }
    }

    #[inline]
    pub fn outer_product(self, rhs: Self) -> SimdMat2 {
        SimdMat2 {
            x_axis: self * rhs.x,
            y_axis: self * rhs.y,
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

impl PointOps for SimdVec2 {
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

impl Add<SimdVec2> for SimdVec2 {
    type Output = Self;
    #[inline]
    fn add(self, rhs: Self) -> Self {
        Self {
            x: self.x.add(rhs.x),
            y: self.y.add(rhs.y),
        }
    }
}

impl AddAssign<SimdVec2> for SimdVec2 {
    #[inline]
    fn add_assign(&mut self, rhs: Self) {
        self.x.add_assign(rhs.x);
        self.y.add_assign(rhs.y);
    }
}

impl Add<SimdReal> for SimdVec2 {
    type Output = Self;
    #[inline]
    fn add(self, rhs: SimdReal) -> Self {
        Self {
            x: self.x.add(rhs),
            y: self.y.add(rhs),
        }
    }
}

impl AddAssign<SimdReal> for SimdVec2 {
    #[inline]
    fn add_assign(&mut self, rhs: SimdReal) {
        self.x.add_assign(rhs);
        self.y.add_assign(rhs);
    }
}

impl Add<SimdVec2> for SimdReal {
    type Output = SimdVec2;
    #[inline]
    fn add(self, rhs: SimdVec2) -> SimdVec2 {
        SimdVec2 {
            x: self.add(rhs.x),
            y: self.add(rhs.y),
        }
    }
}

impl Sub<SimdVec2> for SimdVec2 {
    type Output = Self;
    #[inline]
    fn sub(self, rhs: Self) -> Self {
        Self {
            x: self.x.sub(rhs.x),
            y: self.y.sub(rhs.y),
        }
    }
}

impl SubAssign<SimdVec2> for SimdVec2 {
    #[inline]
    fn sub_assign(&mut self, rhs: SimdVec2) {
        self.x.sub_assign(rhs.x);
        self.y.sub_assign(rhs.y);
    }
}

impl Sub<SimdReal> for SimdVec2 {
    type Output = Self;
    #[inline]
    fn sub(self, rhs: SimdReal) -> Self {
        Self {
            x: self.x.sub(rhs),
            y: self.y.sub(rhs),
        }
    }
}

impl SubAssign<SimdReal> for SimdVec2 {
    #[inline]
    fn sub_assign(&mut self, rhs: SimdReal) {
        self.x.sub_assign(rhs);
        self.y.sub_assign(rhs);
    }
}

impl Sub<SimdVec2> for SimdReal {
    type Output = SimdVec2;
    #[inline]
    fn sub(self, rhs: SimdVec2) -> SimdVec2 {
        SimdVec2 {
            x: self.sub(rhs.x),
            y: self.sub(rhs.y),
        }
    }
}

impl Mul<SimdVec2> for SimdVec2 {
    type Output = Self;
    #[inline]
    fn mul(self, rhs: Self) -> Self {
        Self {
            x: self.x.mul(rhs.x),
            y: self.y.mul(rhs.y),
        }
    }
}

impl MulAssign<SimdVec2> for SimdVec2 {
    #[inline]
    fn mul_assign(&mut self, rhs: Self) {
        self.x.mul_assign(rhs.x);
        self.y.mul_assign(rhs.y);
    }
}

impl<'a> Mul<SimdReal> for &'a SimdVec2 {
    type Output = SimdVec2;
    #[inline]
    fn mul(self, rhs: SimdReal) -> SimdVec2 {
        SimdVec2 {
            x: self.x.mul(rhs),
            y: self.y.mul(rhs),
        }
    }
}

impl Mul<SimdReal> for SimdVec2 {
    type Output = Self;
    #[inline]
    fn mul(self, rhs: SimdReal) -> Self {
        Self {
            x: self.x.mul(rhs),
            y: self.y.mul(rhs),
        }
    }
}

impl MulAssign<SimdReal> for SimdVec2 {
    #[inline]
    fn mul_assign(&mut self, rhs: SimdReal) {
        self.x.mul_assign(rhs);
        self.y.mul_assign(rhs);
    }
}

impl Mul<SimdVec2> for SimdReal {
    type Output = SimdVec2;
    #[inline]
    fn mul(self, rhs: SimdVec2) -> SimdVec2 {
        SimdVec2 {
            x: self.mul(rhs.x),
            y: self.mul(rhs.y),
        }
    }
}

impl Div<SimdVec2> for SimdVec2 {
    type Output = Self;
    #[inline]
    fn div(self, rhs: Self) -> Self {
        Self {
            x: self.x.div(rhs.x),
            y: self.y.div(rhs.y),
        }
    }
}

impl DivAssign<SimdVec2> for SimdVec2 {
    #[inline]
    fn div_assign(&mut self, rhs: Self) {
        self.x.div_assign(rhs.x);
        self.y.div_assign(rhs.y);
    }
}

impl Div<SimdReal> for SimdVec2 {
    type Output = Self;
    #[inline]
    fn div(self, rhs: SimdReal) -> Self {
        Self {
            x: self.x.div(rhs),
            y: self.y.div(rhs),
        }
    }
}

impl DivAssign<SimdReal> for SimdVec2 {
    #[inline]
    fn div_assign(&mut self, rhs: SimdReal) {
        self.x.div_assign(rhs);
        self.y.div_assign(rhs);
    }
}

impl Div<SimdVec2> for SimdReal {
    type Output = SimdVec2;
    #[inline]
    fn div(self, rhs: SimdVec2) -> SimdVec2 {
        SimdVec2 {
            x: self.div(rhs.x),
            y: self.div(rhs.y),
        }
    }
}

impl Rem<SimdVec2> for SimdVec2 {
    type Output = Self;
    #[inline]
    fn rem(self, rhs: Self) -> Self {
        Self {
            x: self.x.rem(rhs.x),
            y: self.y.rem(rhs.y),
        }
    }
}

impl RemAssign<SimdVec2> for SimdVec2 {
    #[inline]
    fn rem_assign(&mut self, rhs: Self) {
        self.x.rem_assign(rhs.x);
        self.y.rem_assign(rhs.y);
    }
}

impl Rem<SimdReal> for SimdVec2 {
    type Output = Self;
    #[inline]
    fn rem(self, rhs: SimdReal) -> Self {
        Self {
            x: self.x.rem(rhs),
            y: self.y.rem(rhs),
        }
    }
}

impl RemAssign<SimdReal> for SimdVec2 {
    #[inline]
    fn rem_assign(&mut self, rhs: SimdReal) {
        self.x.rem_assign(rhs);
        self.y.rem_assign(rhs);
    }
}

impl Rem<SimdVec2> for SimdReal {
    type Output = SimdVec2;
    #[inline]
    fn rem(self, rhs: SimdVec2) -> SimdVec2 {
        SimdVec2 {
            x: self.rem(rhs.x),
            y: self.rem(rhs.y),
        }
    }
}

impl Neg for SimdVec2 {
    type Output = Self;
    #[inline]
    fn neg(self) -> Self::Output {
        Self {
            x: -self.x,
            y: -self.y,
        }
    }
}

#[derive(Clone, Copy, Default, Debug, PartialEq)]
pub struct SimdIso2 {
    pub translation: SimdVec2,
    pub rotation: SimdMat2,
}

impl SimdIso2 {
    // /// An identity isometry.
    // pub const IDENTITY: Self = Self {
    //     translation: SimdVec2::zeros(),
    //     rotation: SimdQuat::IDENTITY,
    // };

    #[inline]
    pub fn new(translation: Vec2, rotation: Real) -> Self {
        Self {
            translation: SimdVec2::splat(translation),
            rotation: SimdMat2::from_angle(rotation),
        }
    }

    #[inline]
    pub fn identity() -> Self {
        Self {
            translation: SimdVec2::zeros(),
            rotation: SimdMat2::identity(),
        }
    }

    #[inline]
    pub fn splat(iso: Iso2) -> Self {
        Self {
            translation: SimdVec2::splat(iso.translation),
            rotation: SimdMat2::splat(iso.rotation),
        }
    }

    #[inline]
    pub fn extract(&self, lane: usize) -> Iso2 {
        Iso2 {
            translation: self.translation.extract(lane),
            rotation: self.rotation.extract(lane),
        }
    }

    #[inline]
    pub fn from_translation(translation: Vec2) -> Self {
        Self {
            translation: SimdVec2::splat(translation),
            rotation: SimdMat2::identity(),
        }
    }

    #[inline]
    pub fn from_rotation(rotation: Mat2) -> Self {
        Self {
            translation: SimdVec2::zeros(),
            rotation: SimdMat2::splat(rotation),
        }
    }

    #[inline]
    pub fn inverse(self) -> Self {
        let inv_rot = self.rotation.transpose();
        Self {
            translation: inv_rot * -self.translation,
            rotation: inv_rot,
        }
    }

    #[inline]
    pub fn transform_point(self, point: SimdVec2) -> SimdVec2 {
        self.translation + self.rotation * point
    }

    #[inline]
    pub fn inverse_transform_point(self, point: SimdVec2) -> SimdVec2 {
        self.rotation.transpose() * (point - self.translation)
    }
}

impl From<[Iso2; 4]> for SimdIso2 {
    fn from(value: [Iso2; 4]) -> Self {
        Self {
            rotation: value.map(|iso| iso.rotation).into(),
            translation: value.map(|tra| tra.translation).into(),
        }
    }
}

impl Mul<SimdIso2> for SimdIso2 {
    type Output = SimdIso2;
    #[inline]
    fn mul(self, rhs: SimdIso2) -> Self::Output {
        Self {
            rotation: self.rotation * rhs.rotation,
            translation: self.translation + self.rotation * rhs.translation,
        }
    }
}

impl Mul<SimdVec2> for SdpMatrix2<SimdReal> {
    type Output = SimdVec2;

    #[inline]
    fn mul(self, rhs: SimdVec2) -> Self::Output {
        let x = self.m11 * rhs.x + self.m12 * rhs.y;
        let y = self.m12 * rhs.x + self.m22 * rhs.y;
        SimdVec2::new(x, y)
    }
}
