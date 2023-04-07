use core::ops::IndexMut;
#[cfg(feature = "std")]
use na::{DMatrix, DVector, Scalar};

/// Abstraction over a 1D array.
pub trait Array1<T>: IndexMut<usize, Output = T> {
    /// The number of heights on this storage.
    fn len(&self) -> usize;

    /// Is this array empty?
    #[inline]
    fn is_empty(&self) -> bool {
        self.len() == 0
    }

    // NOTE: we don’t name it just `get` to avoid clashes with pre-existing `get` methods.
    /// Gets the i-th element of this array, if it exists.
    #[inline]
    fn get_at(&self, i: usize) -> Option<&T> {
        if i < self.len() {
            Some(&self[i])
        } else {
            None
        }
    }
}

#[cfg(feature = "std")]
impl<T> Array1<T> for Vec<T> {
    #[inline]
    fn len(&self) -> usize {
        self.len()
    }
}

#[cfg(feature = "std")]
impl<T> Array1<T> for DVector<T> {
    #[inline]
    fn len(&self) -> usize {
        self.len()
    }
}

/// Abstraction over a 2D array.
pub trait Array2 {
    /// The type of heights.
    type Item;
    /// The number of rows of the heights grid.
    fn nrows(&self) -> usize;
    /// The number of columns of the heights grid.
    fn ncols(&self) -> usize;
    /// Gets the height on the `(i, j)`-th cell of the height grid.
    fn get(&self, i: usize, j: usize) -> Self::Item;
    /// Sets the height on the `(i, j)`-th cell of the height grid.
    fn set(&mut self, i: usize, j: usize, val: Self::Item);
}

#[cfg(feature = "std")]
impl<T: Scalar> Array2 for DMatrix<T> {
    type Item = T;

    #[inline]
    fn nrows(&self) -> usize {
        self.nrows()
    }

    #[inline]
    fn ncols(&self) -> usize {
        self.ncols()
    }

    #[inline]
    fn get(&self, i: usize, j: usize) -> Self::Item {
        self[(i, j)].clone()
    }

    #[inline]
    fn set(&mut self, i: usize, j: usize, val: Self::Item) {
        self[(i, j)] = val
    }
}

#[cfg_attr(feature = "serde-serialize", derive(Serialize, Deserialize))]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Deserialize, rkyv::Serialize),
    archive(check_bytes)
)]
#[derive(Copy, Clone, Debug, PartialEq, Eq, Default)]
/// Default data storage based on `Vec`, `DVector`, and `DMatrix`.
pub struct DefaultStorage;
#[derive(Copy, Clone, Debug, PartialEq, Eq, Default)]
#[cfg(feature = "cuda")]
/// Data storage residing in CUDA memory.
pub struct CudaStorage;
#[derive(Copy, Clone, Debug, PartialEq, Eq, Default)]
#[cfg(feature = "cuda")]
/// Data storage for accessing data from CUDA kernels.
pub struct CudaStoragePtr;
