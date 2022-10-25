use crate::math::{Point, Real, Vector};
use crate::partitioning::QBVHStorage;
use crate::shape::{TopoFace, TopoHalfEdge, TopoVertex};
use crate::utils::{Array1, DefaultStorage};

#[cfg(all(feature = "std", feature = "cuda"))]
use crate::utils::CudaArray1;
#[cfg(feature = "cuda")]
use crate::utils::{CudaArrayPointer1, CudaStorage, CudaStoragePtr};

pub trait TriMeshStorage {
    type QBVHStorage: QBVHStorage<u32>;
    type ArrayTopoVertex: Array1<TopoVertex>;
    type ArrayTopoFace: Array1<TopoFace>;
    type ArrayTopoHalfEdge: Array1<TopoHalfEdge>;
    type ArrayU32: Array1<u32>;
    type ArrayUsize: Array1<usize>;
    type ArrayVector: Array1<Vector<Real>>;
    type ArrayPoint: Array1<Point<Real>>;
    type ArrayIdx: Array1<[u32; 3]>;
    type ArrayVectorTriple: Array1<[Vector<Real>; 3]>;
}

#[cfg(feature = "std")]
impl TriMeshStorage for DefaultStorage {
    type QBVHStorage = Self;
    type ArrayTopoVertex = Vec<TopoVertex>;
    type ArrayTopoFace = Vec<TopoFace>;
    type ArrayTopoHalfEdge = Vec<TopoHalfEdge>;
    type ArrayU32 = Vec<u32>;
    type ArrayUsize = Vec<usize>;
    type ArrayVector = Vec<Vector<Real>>;
    type ArrayPoint = Vec<Point<Real>>;
    type ArrayIdx = Vec<[u32; 3]>;
    type ArrayVectorTriple = Vec<[Vector<Real>; 3]>;
}

#[cfg(all(feature = "std", feature = "cuda"))]
impl TriMeshStorage for CudaStorage {
    type QBVHStorage = Self;
    type ArrayTopoVertex = CudaArray1<TopoVertex>;
    type ArrayTopoFace = CudaArray1<TopoFace>;
    type ArrayTopoHalfEdge = CudaArray1<TopoHalfEdge>;
    type ArrayU32 = CudaArray1<u32>;
    type ArrayUsize = CudaArray1<usize>;
    type ArrayVector = CudaArray1<Vector<Real>>;
    type ArrayPoint = CudaArray1<Point<Real>>;
    type ArrayIdx = CudaArray1<[u32; 3]>;
    type ArrayVectorTriple = CudaArray1<[Vector<Real>; 3]>;
}

#[cfg(feature = "cuda")]
impl TriMeshStorage for CudaStoragePtr {
    type QBVHStorage = Self;
    type ArrayTopoVertex = CudaArrayPointer1<TopoVertex>;
    type ArrayTopoFace = CudaArrayPointer1<TopoFace>;
    type ArrayTopoHalfEdge = CudaArrayPointer1<TopoHalfEdge>;
    type ArrayU32 = CudaArrayPointer1<u32>;
    type ArrayUsize = CudaArrayPointer1<usize>;
    type ArrayVector = CudaArrayPointer1<Vector<Real>>;
    type ArrayPoint = CudaArrayPointer1<Point<Real>>;
    type ArrayIdx = CudaArrayPointer1<[u32; 3]>;
    type ArrayVectorTriple = CudaArrayPointer1<[Vector<Real>; 3]>;
}
