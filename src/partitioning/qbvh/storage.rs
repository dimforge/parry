use crate::partitioning::{QBVHNode, QBVHProxy};
use crate::utils::{Array1, DefaultStorage};

#[cfg(all(feature = "std", feature = "cuda"))]
use crate::utils::CudaArray1;
#[cfg(feature = "cuda")]
use crate::utils::{CudaArrayPointer1, CudaStorage, CudaStoragePtr};

/// Trait describing all the types needed for storing a QBVH’s data.
pub trait QBVHStorage<LeafData> {
    /// Type of the array containing the QBVH nodes.
    type Nodes: Array1<QBVHNode>;
    /// Type of an array containing u32.
    type ArrayU32: Array1<u32>;
    /// Type of the array containing the QBVH leaves.
    type ArrayProxies: Array1<QBVHProxy<LeafData>>;
}

#[cfg(feature = "std")]
impl<LeafData> QBVHStorage<LeafData> for DefaultStorage {
    type Nodes = Vec<QBVHNode>;
    type ArrayU32 = Vec<u32>;
    type ArrayProxies = Vec<QBVHProxy<LeafData>>;
}

#[cfg(all(feature = "std", feature = "cuda"))]
impl<LeafData: cust_core::DeviceCopy> QBVHStorage<LeafData> for CudaStorage {
    type Nodes = CudaArray1<QBVHNode>;
    type ArrayU32 = CudaArray1<u32>;
    type ArrayProxies = CudaArray1<QBVHProxy<LeafData>>;
}

#[cfg(feature = "cuda")]
impl<LeafData: cust_core::DeviceCopy> QBVHStorage<LeafData> for CudaStoragePtr {
    type Nodes = CudaArrayPointer1<QBVHNode>;
    type ArrayU32 = CudaArrayPointer1<u32>;
    type ArrayProxies = CudaArrayPointer1<QBVHProxy<LeafData>>;
}
