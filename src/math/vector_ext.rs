//! Vector extension trait for glam types.

use crate::math::{IVector, Int, Matrix, Real, Vector};

#[cfg(not(feature = "std"))]
use simba::scalar::ComplexField;

/// Extension trait for glam vector types to provide additional functionality.
pub trait VectorExt: Sized + Copy {
    /// Creates a vector with the i-th component set to `val` and all others to zero.
    fn ith(i: usize, val: Real) -> Self;

    /// Computes the angle between two vectors in radians.
    fn angle(self, other: Self) -> Real;

    /// Computes the kronecker product between two vectors.
    fn kronecker(self, other: Self) -> Matrix;

    /// Gets the i-th component of the vector by index.
    ///
    /// This is a SPIR-V compatible replacement for bracket indexing (`v[i]`),
    /// which is not supported on glam types when targeting SPIR-V.
    fn vget(&self, i: usize) -> Real;

    /// Sets the i-th component of the vector by index.
    ///
    /// This is a SPIR-V compatible replacement for bracket indexing (`v[i] = val`),
    /// which is not supported on glam types when targeting SPIR-V.
    fn vset(&mut self, i: usize, val: Real);
}

impl VectorExt for Vector {
    #[inline]
    #[cfg(feature = "dim2")]
    fn ith(i: usize, val: Real) -> Self {
        match i {
            0 => Self::new(val, 0.0),
            1 => Self::new(0.0, val),
            _ => Self::ZERO,
        }
    }

    #[inline]
    #[cfg(feature = "dim3")]
    fn ith(i: usize, val: Real) -> Self {
        match i {
            0 => Self::new(val, 0.0, 0.0),
            1 => Self::new(0.0, val, 0.0),
            2 => Self::new(0.0, 0.0, val),
            _ => Self::ZERO,
        }
    }

    #[inline]
    fn angle(self, other: Self) -> Real {
        let dot = self.dot(other);
        let len_self = self.length();
        let len_other = other.length();
        if len_self > 0.0 && len_other > 0.0 {
            (dot / (len_self * len_other)).clamp(-1.0, 1.0).acos()
        } else {
            0.0
        }
    }

    #[inline]
    #[cfg(feature = "dim2")]
    fn kronecker(self, other: Self) -> Matrix {
        Matrix::from_cols(self * other.x, self * other.y)
    }

    #[inline]
    #[cfg(feature = "dim3")]
    fn kronecker(self, other: Self) -> Matrix {
        Matrix::from_cols(self * other.x, self * other.y, self * other.z)
    }

    #[inline]
    fn vget(&self, i: usize) -> Real {
        #[cfg(target_arch = "spirv")]
        return match i {
            0 => self.x,
            1 => self.y,
            #[cfg(feature = "dim3")]
            2 => self.z,
            _ => panic!("Index out of bounds for 2D vector: {i}"),
        };
        #[cfg(not(target_arch = "spirv"))]
        {
            self[i]
        }
    }

    #[inline]
    fn vset(&mut self, i: usize, val: Real) {
        #[cfg(target_arch = "spirv")]
        match i {
            0 => self.x = val,
            1 => self.y = val,
            #[cfg(feature = "dim3")]
            2 => self.z = val,
            _ => panic!("Index out of bounds for 2D vector: {i}"),
        };
        #[cfg(not(target_arch = "spirv"))]
        {
            self[i] = val;
        }
    }
}

/// Extension trait for integer vector types with index-based component access.
///
/// Provides SPIR-V compatible replacement for bracket indexing on glam integer vectors.
pub trait IVectorExt: Sized + Copy {
    /// Gets the i-th component of the integer vector by index.
    fn ivget(&self, i: usize) -> Int;

    /// Sets the i-th component of the integer vector by index.
    fn ivset(&mut self, i: usize, val: Int);
}

impl IVectorExt for IVector {
    #[inline]
    fn ivget(&self, i: usize) -> Int {
        #[cfg(target_arch = "spirv")]
        return match i {
            0 => self.x,
            1 => self.y,
            #[cfg(feature = "dim3")]
            2 => self.z,
            _ => panic!("Index out of bounds for 2D integer vector: {i}"),
        };
        #[cfg(not(target_arch = "spirv"))]
        {
            self[i]
        }
    }

    #[inline]
    fn ivset(&mut self, i: usize, val: Int) {
        #[cfg(target_arch = "spirv")]
        match i {
            0 => self.x = val,
            1 => self.y = val,
            #[cfg(feature = "dim3")]
            2 => self.z = val,
            _ => panic!("Index out of bounds for 2D integer vector: {i}"),
        }
        #[cfg(not(target_arch = "spirv"))]
        {
            self[i] = val;
        }
    }
}
