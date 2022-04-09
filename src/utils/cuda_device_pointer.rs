#[cfg(not(target_os = "cuda"))]
pub use cust::memory::DevicePointer;

#[cfg(target_os = "cuda")]
pub use self::cuda_utils::*;

#[cfg(target_os = "cuda")]
pub mod cuda_utils {
    #[repr(transparent)]
    #[derive(PartialEq)]
    pub struct DevicePointer<T>(*mut T);

    impl<T> Clone for DevicePointer<T> {
        fn clone(&self) -> Self {
            Self(self.0)
        }
    }

    impl<T> Copy for DevicePointer<T> {}
    unsafe impl<T> cust_core::DeviceCopy for DevicePointer<T> {}

    impl<T> DevicePointer<T> {
        pub fn null() -> Self {
            Self(core::ptr::null_mut())
        }
        pub fn as_ptr(self) -> *const T {
            self.0
        }
        pub fn as_mut_ptr(&self) -> *mut T {
            self.0
        }
    }
}
