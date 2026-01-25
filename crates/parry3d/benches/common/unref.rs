use parry3d::math::{Pose3, Vec3};

pub trait Unref<T> {
    fn unref(a: Self) -> T;
}

impl<'a> Unref<f32> for &'a f32 {
    #[inline(always)]
    fn unref(a: &f32) -> f32 {
        *a
    }
}

impl<'a> Unref<f64> for &'a f64 {
    #[inline(always)]
    fn unref(a: &f64) -> f64 {
        *a
    }
}

impl<'a> Unref<bool> for &'a bool {
    #[inline(always)]
    fn unref(a: &bool) -> bool {
        *a
    }
}

impl<'a> Unref<Vec3> for &'a Vec3 {
    #[inline(always)]
    fn unref(a: &Vec3) -> Vec3 {
        *a
    }
}

impl<'a> Unref<Pose3> for &'a Pose3 {
    #[inline(always)]
    fn unref(a: &Pose3) -> Pose3 {
        *a
    }
}

impl<'a, T> Unref<&'a T> for &'a T {
    #[inline(always)]
    fn unref(a: &'a T) -> &'a T {
        a
    }
}

#[inline(always)]
pub fn unref<T: Unref<O>, O>(val: T) -> O {
    Unref::unref(val)
}
