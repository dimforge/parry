#[cfg(target_os = "cuda")]
use crate::shape::HeightFieldStorage;
use crate::utils::DevicePointer;

#[cfg(feature = "std")]
use cust::{error::CudaResult, memory::DeviceBuffer};
use cust_core::DeviceCopy;

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

#[repr(C)]
#[derive(Copy, Clone, cust_core::DeviceCopy)]
/// A pointer to a 2D CUDA array.
pub struct CudaArrayPointer2<T: ?Sized + DeviceCopy> {
    data: DevicePointer<T>,
    nrows: usize,
    ncols: usize,
}

#[cfg(all(feature = "dim3", target_os = "cuda"))]
impl<T: ?Sized + DeviceCopy> HeightFieldStorage for CudaArrayPointer2<T> {
    type Item = T;

    fn nrows(&self) -> usize {
        self.nrows
    }

    fn ncols(&self) -> usize {
        self.ncols
    }

    fn get(&self, i: usize, j: usize) -> Self::Item {
        let linear_index = i + j * self.nrows;
        assert!(linear_index < self.nrows * self.ncols);
        unsafe { *self.data.as_ptr().add(linear_index) }
    }

    fn set(&mut self, i: usize, j: usize, val: Self::Item) {
        let linear_index = i + j * self.nrows;
        assert!(linear_index < self.nrows * self.ncols);
        unsafe {
            *self.data.as_mut_ptr().add(linear_index) = val;
        }
    }
}

/*
 *
 * 1D array.
 *
 */
#[cfg(feature = "std")]
/// A 1D array residing on GPU memory.
pub struct CudaArray1<T: ?Sized + DeviceCopy> {
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

#[repr(C)]
#[derive(Copy, Clone, cust_core::DeviceCopy)]
/// A pointer to a 2D CUDA array.
pub struct CudaArrayPointer1<T: ?Sized + DeviceCopy> {
    data: DevicePointer<T>,
    len: usize,
}

#[cfg(all(feature = "dim2", target_os = "cuda"))]
impl<T: ?Sized + DeviceCopy> HeightFieldStorage for CudaArrayPointer1<T> {
    type Item = T;

    fn len(&self) -> usize {
        self.len
    }

    fn get(&self, i: usize) -> Self::Item {
        assert!(i < self.len);
        unsafe { *self.data.as_ptr().add(i) }
    }

    fn set(&mut self, i: usize, val: Self::Item) {
        assert!(i < self.len);
        unsafe {
            *self.data.as_mut_ptr().add(i) = val;
        }
    }
}
