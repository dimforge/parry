#[cfg(target_os = "cuda")]
use crate::shape::HeightFieldStorage;
use crate::utils::{Array1, Array2, DevicePointer};

#[cfg(feature = "std")]
use cust::{error::CudaResult, memory::DeviceBuffer};
use cust_core::DeviceCopy;
use std::ops::{Index, IndexMut};

/*
 *
 * 2D array.
 *
*/
#[cfg(feature = "std")]
/// A 2D array residing on GPU memory.
pub struct CudaArray2<T: ?Sized + DeviceCopy> {
    data: DeviceBuffer<T>,
    nrows: usize,
    ncols: usize,
}

#[cfg(feature = "std")]
impl<T: ?Sized + DeviceCopy> CudaArray2<T> {
    /// Initialize a 2D cuda array on the GPU.
    pub fn new(data: &[T], nrows: usize, ncols: usize) -> CudaResult<Self> {
        assert_eq!(
            data.len(),
            nrows * ncols,
            "The data length much match the array length."
        );
        DeviceBuffer::from_slice(data).map(|data| Self { data, nrows, ncols })
    }

    /// Initialize, using a matrix, a 2D cuda array on the GPU.
    pub fn from_matrix(mat: &na::DMatrix<T>) -> CudaResult<Self> {
        Self::new(mat.as_slice(), mat.nrows(), mat.ncols())
    }

    /// Gets the device pointer to the CUDA memory.
    pub fn as_device_ptr(&self) -> CudaArrayPointer2<T> {
        CudaArrayPointer2 {
            data: self.data.as_device_ptr(),
            nrows: self.nrows,
            ncols: self.ncols,
        }
    }
}

#[cfg(feature = "std")]
impl<T: ?Sized + DeviceCopy> Array2 for CudaArray2<T> {
    type Item = T;

    fn nrows(&self) -> usize {
        panic!("Cuda arrays cannot be read directly.");
    }

    fn ncols(&self) -> usize {
        panic!("Cuda arrays cannot be read directly.");
    }

    fn get(&self, i: usize, j: usize) -> Self::Item {
        panic!("Cuda arrays cannot be read directly.");
    }

    fn set(&mut self, i: usize, j: usize, val: Self::Item) {
        panic!("Cuda arrays cannot be read directly.");
    }
}

#[repr(C)]
#[derive(Copy, Clone, cust_core::DeviceCopy)]
/// A pointer to a 2D CUDA array.
pub struct CudaArrayPointer2<T: ?Sized + DeviceCopy> {
    data: DevicePointer<T>,
    nrows: usize,
    ncols: usize,
}

#[cfg(all(feature = "dim3", target_os = "cuda"))]
impl<T: ?Sized + DeviceCopy> Array2 for CudaArrayPointer2<T> {
    type Item = T;

    #[inline]
    fn nrows(&self) -> usize {
        self.nrows
    }

    #[inline]
    fn ncols(&self) -> usize {
        self.ncols
    }

    #[inline]
    fn get(&self, i: usize, j: usize) -> Self::Item {
        let linear_index = i + j * self.nrows;
        assert!(linear_index < self.nrows * self.ncols);
        unsafe { *self.data.as_ptr().add(linear_index) }
    }

    #[inline]
    fn set(&mut self, i: usize, j: usize, val: Self::Item) {
        let linear_index = i + j * self.nrows;
        assert!(linear_index < self.nrows * self.ncols);
        unsafe {
            *self.data.as_mut_ptr().add(linear_index) = val;
        }
    }
}

#[cfg(all(feature = "dim3", not(target_os = "cuda")))]
impl<T: ?Sized + DeviceCopy> Array2 for CudaArrayPointer2<T> {
    type Item = T;

    fn nrows(&self) -> usize {
        panic!("Cuda pointers can only be read from of a cuda kernel.");
    }

    fn ncols(&self) -> usize {
        panic!("Cuda pointers can only be read from of a cuda kernel.");
    }

    fn get(&self, i: usize, j: usize) -> Self::Item {
        panic!("Cuda pointers can only be read from of a cuda kernel.");
    }

    fn set(&mut self, i: usize, j: usize, val: Self::Item) {
        panic!("Cuda pointers can only be read from of a cuda kernel.");
    }
}

/*
 *
 * 1D array.
 *
 */
#[cfg(feature = "std")]
/// A 1D array residing on GPU memory.
pub struct CudaArray1<T: DeviceCopy> {
    data: DeviceBuffer<T>,
}

#[cfg(feature = "std")]
impl<T: ?Sized + DeviceCopy> CudaArray1<T> {
    /// Initialize a 1D cuda array on the GPU.
    pub fn new(data: &[T]) -> CudaResult<Self> {
        DeviceBuffer::from_slice(data).map(|data| Self { data })
    }

    /// Initialize a 1D cuda array on the GPU using a dynamically-sized vector.
    pub fn from_vector(vect: &na::DVector<T>) -> CudaResult<Self> {
        Self::new(vect.as_slice())
    }

    /// Gets the device pointer to the CUDA memory.
    pub fn as_device_ptr(&self) -> CudaArrayPointer1<T> {
        CudaArrayPointer1 {
            data: self.data.as_device_ptr(),
            len: self.data.len(),
        }
    }
}

#[cfg(feature = "std")]
impl<T: ?Sized + DeviceCopy> Array1<T> for CudaArray1<T> {
    fn len(&self) -> usize {
        panic!("Cuda arrays cannot be read directly.");
    }
}

#[cfg(feature = "std")]
impl<T: DeviceCopy> Index<usize> for CudaArray1<T> {
    type Output = T;

    #[inline]
    fn index(&self, _: usize) -> &T {
        panic!("Cuda arrays cannot be read directly.");
    }
}

#[cfg(feature = "std")]
impl<T: DeviceCopy> IndexMut<usize> for CudaArray1<T> {
    #[inline]
    fn index_mut(&mut self, _: usize) -> &mut T {
        panic!("Cuda arrays cannot be read directly.");
    }
}

#[repr(C)]
#[derive(Copy, Clone, cust_core::DeviceCopy)]
/// A pointer to a 2D CUDA array.
pub struct CudaArrayPointer1<T: ?Sized + DeviceCopy> {
    data: DevicePointer<T>,
    len: usize,
}

#[cfg(target_os = "cuda")]
impl<T: ?Sized + DeviceCopy> CudaArrayPointer1<T> {
    #[inline]
    pub fn len(&self) -> usize {
        self.len
    }

    #[inline]
    pub fn get(&self, i: usize) -> T {
        assert!(i < self.len);
        unsafe { *self.data.as_ptr().add(i) }
    }

    #[inline]
    pub fn set(&mut self, i: usize, val: T) {
        assert!(i < self.len);
        unsafe {
            *self.data.as_mut_ptr().add(i) = val;
        }
    }
}

#[cfg(not(target_os = "cuda"))]
impl<T: ?Sized + DeviceCopy> Array1<T> for CudaArrayPointer1<T> {
    fn len(&self) -> usize {
        panic!("Cuda pointers can only be read from of a cuda kernel.");
    }
}

#[cfg(target_os = "cuda")]
impl<T: ?Sized + DeviceCopy> Array1<T> for CudaArrayPointer1<T> {
    #[inline]
    fn len(&self) -> usize {
        self.len
    }
}

#[cfg(not(target_os = "cuda"))]
impl<T: ?Sized + DeviceCopy> Index<usize> for CudaArrayPointer1<T> {
    type Output = T;

    #[inline]
    fn index(&self, i: usize) -> &T {
        panic!("Cuda pointers can only be read from of a cuda kernel.");
    }
}

#[cfg(not(target_os = "cuda"))]
impl<T: ?Sized + DeviceCopy> IndexMut<usize> for CudaArrayPointer1<T> {
    #[inline]
    fn index_mut(&mut self, i: usize) -> &mut T {
        panic!("Cuda pointers can only be read from of a cuda kernel.");
    }
}

#[cfg(target_os = "cuda")]
impl<T: DeviceCopy> Index<usize> for CudaArrayPointer1<T> {
    type Output = T;

    #[inline]
    fn index(&self, i: usize) -> &T {
        assert!(i < self.len);
        unsafe { &*self.data.as_ptr().add(i) }
    }
}

#[cfg(target_os = "cuda")]
impl<T: DeviceCopy> IndexMut<usize> for CudaArrayPointer1<T> {
    #[inline]
    fn index_mut(&mut self, i: usize) -> &mut T {
        assert!(i < self.len);
        unsafe { &mut *self.data.as_mut_ptr().add(i) }
    }
}
