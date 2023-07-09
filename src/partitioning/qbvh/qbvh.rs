use crate::bounding_volume::{Aabb, SimdAabb};
use crate::math::{Real, Vector};
use crate::partitioning::qbvh::storage::QbvhStorage;
use crate::utils::DefaultStorage;
use bitflags::bitflags;

use na::SimdValue;

#[cfg(feature = "rkyv")]
use rkyv::{bytecheck, CheckBytes};
#[cfg(all(feature = "std", feature = "cuda"))]
use {crate::utils::CudaArray1, cust::error::CudaResult};
#[cfg(feature = "cuda")]
use {
    crate::utils::{CudaStorage, CudaStoragePtr},
    cust_core::DeviceCopy,
};

/// A data to which an index is associated.
pub trait IndexedData: Copy {
    /// Creates a new default instance of `Self`.
    fn default() -> Self;
    /// Gets the index associated to `self`.
    fn index(&self) -> usize;
}

impl IndexedData for usize {
    fn default() -> Self {
        // NOTE: we use u32::MAX for compatibility
        // between 32 and 64 bit platforms.
        u32::MAX as usize
    }
    fn index(&self) -> usize {
        *self
    }
}

impl IndexedData for u32 {
    fn default() -> Self {
        u32::MAX
    }
    fn index(&self) -> usize {
        *self as usize
    }
}

impl IndexedData for u64 {
    fn default() -> Self {
        u64::MAX as u64
    }
    fn index(&self) -> usize {
        *self as usize
    }
}

/// The index of an internal SIMD node of a Qbvh.
pub type SimdNodeIndex = u32;

/// The index of a node part of a Qbvh.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
#[cfg_attr(feature = "serde-serialize", derive(Serialize, Deserialize))]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Deserialize, rkyv::Serialize, CheckBytes),
    archive(as = "Self")
)]
#[cfg_attr(feature = "cuda", derive(cust_core::DeviceCopy))]
/// The index of one specific node of a Qbvh.
pub struct NodeIndex {
    /// The index of the SIMD node containing the addressed node.
    pub index: SimdNodeIndex, // Index of the addressed node in the `nodes` array.
    /// The SIMD lane the addressed node is associated to.
    pub lane: u8, // SIMD lane of the addressed node.
}

impl NodeIndex {
    pub(super) fn new(index: u32, lane: u8) -> Self {
        Self { index, lane }
    }

    pub(super) fn invalid() -> Self {
        Self {
            index: u32::MAX,
            lane: 0,
        }
    }

    pub(super) fn is_invalid(&self) -> bool {
        self.index == u32::MAX
    }
}

bitflags! {
    #[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
    #[cfg_attr(
        feature = "rkyv",
        derive(rkyv::Archive, rkyv::Deserialize, rkyv::Serialize),
        archive(as = "Self")
    )]
    #[cfg_attr(feature = "cuda", derive(cust_core::DeviceCopy))]
    #[derive(Default)]
    /// The status of a QBVH node.
    pub struct QbvhNodeFlags: u8 {
        /// If this bit is set, the node is a leaf.
        const LEAF = 0b0001;
        /// If this bit is set, this node was recently changed.
        const CHANGED = 0b0010;
        /// Does this node need an update?
        const DIRTY = 0b0100;
    }
}

/// A SIMD node of an SIMD Qbvh.
///
/// This groups four nodes of the Qbvh.
#[derive(Copy, Clone, Debug)]
#[cfg_attr(feature = "serde-serialize", derive(Serialize, Deserialize))]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Deserialize, rkyv::Serialize),
    archive(check_bytes)
)]
#[cfg_attr(feature = "cuda", derive(cust_core::DeviceCopy))]
pub struct QbvhNode {
    /// The Aabbs of the qbvh nodes represented by this node.
    pub simd_aabb: SimdAabb,
    /// Index of the nodes of the 4 nodes represented by `self`.
    /// If this is a leaf, it contains the proxy ids instead.
    pub children: [u32; 4],
    /// The index of the node parent to the 4 nodes represented by `self`.
    pub parent: NodeIndex,
    /// Status flags for this node.
    pub flags: QbvhNodeFlags,
}

impl QbvhNode {
    #[inline]
    /// Is this node a leaf?
    pub fn is_leaf(&self) -> bool {
        self.flags.contains(QbvhNodeFlags::LEAF)
    }

    #[inline]
    /// Does the AABB of this node needs to be updated?
    pub fn is_dirty(&self) -> bool {
        self.flags.contains(QbvhNodeFlags::DIRTY)
    }

    #[inline]
    /// Sets if the AABB of this node needs to be updated.
    pub fn set_dirty(&mut self, dirty: bool) {
        self.flags.set(QbvhNodeFlags::DIRTY, dirty);
    }

    #[inline]
    /// Was the AABB of this node changed since the last rebalancing?
    pub fn is_changed(&self) -> bool {
        self.flags.contains(QbvhNodeFlags::CHANGED)
    }

    #[inline]
    /// Sets if the AABB of this node changed since the last rebalancing.
    pub fn set_changed(&mut self, changed: bool) {
        self.flags.set(QbvhNodeFlags::CHANGED, changed);
    }

    #[inline]
    /// An empty internal node.
    pub fn empty() -> Self {
        Self {
            simd_aabb: SimdAabb::new_invalid(),
            children: [u32::MAX; 4],
            parent: NodeIndex::invalid(),
            flags: QbvhNodeFlags::default(),
        }
    }

    #[inline]
    /// An empty leaf.
    pub fn empty_leaf_with_parent(parent: NodeIndex) -> Self {
        Self {
            simd_aabb: SimdAabb::new_invalid(),
            children: [u32::MAX; 4],
            parent,
            flags: QbvhNodeFlags::LEAF,
        }
    }
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
#[cfg_attr(feature = "serde-serialize", derive(Serialize, Deserialize))]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Deserialize, rkyv::Serialize),
    archive(check_bytes)
)]
#[cfg_attr(feature = "cuda", derive(cust_core::DeviceCopy))]
/// Combination of a leaf data and its associated node’s index.
pub struct QbvhProxy<LeafData> {
    /// Index of the leaf node the leaf data is associated to.
    pub node: NodeIndex,
    /// The data contained in this node.
    pub data: LeafData, // The collider data. TODO: only set the collider generation here?
}

impl<LeafData> QbvhProxy<LeafData> {
    pub(super) fn invalid() -> Self
    where
        LeafData: IndexedData,
    {
        Self {
            node: NodeIndex::invalid(),
            data: LeafData::default(),
        }
    }

    pub(super) fn detached(data: LeafData) -> Self {
        Self {
            node: NodeIndex::invalid(),
            data,
        }
    }

    pub(super) fn is_detached(&self) -> bool {
        self.node.is_invalid()
    }
}

/// A quaternary bounding-volume-hierarchy.
///
/// This is a bounding-volume-hierarchy where each node has either four children or none.
#[cfg_attr(feature = "serde-serialize", derive(Serialize, Deserialize))]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Deserialize, rkyv::Serialize),
    archive(check_bytes)
)]
#[repr(C)] // Needed for Cuda.
#[derive(Debug)]
pub struct GenericQbvh<LeafData, Storage: QbvhStorage<LeafData>> {
    pub(super) root_aabb: Aabb,
    pub(super) nodes: Storage::Nodes,
    pub(super) dirty_nodes: Storage::ArrayU32,
    pub(super) free_list: Storage::ArrayU32,
    pub(super) proxies: Storage::ArrayProxies,
}

/// A quaternary bounding-volume-hierarchy.
///
/// This is a bounding-volume-hierarchy where each node has either four children or none.
pub type Qbvh<LeafData> = GenericQbvh<LeafData, DefaultStorage>;
#[cfg(feature = "cuda")]
/// A Qbvh stored in CUDA memory.
pub type CudaQbvh<LeafData> = GenericQbvh<LeafData, CudaStorage>;
#[cfg(feature = "cuda")]
/// A Qbvh accessible from CUDA kernels.
pub type CudaQbvhPtr<LeafData> = GenericQbvh<LeafData, CudaStoragePtr>;

impl<LeafData, Storage> Clone for GenericQbvh<LeafData, Storage>
where
    Storage: QbvhStorage<LeafData>,
    Storage::Nodes: Clone,
    Storage::ArrayU32: Clone,
    Storage::ArrayProxies: Clone,
{
    fn clone(&self) -> Self {
        Self {
            root_aabb: self.root_aabb,
            nodes: self.nodes.clone(),
            dirty_nodes: self.dirty_nodes.clone(),
            free_list: self.free_list.clone(),
            proxies: self.proxies.clone(),
        }
    }
}

impl<LeafData, Storage> Copy for GenericQbvh<LeafData, Storage>
where
    Storage: QbvhStorage<LeafData>,
    Storage::Nodes: Copy,
    Storage::ArrayU32: Copy,
    Storage::ArrayProxies: Copy,
{
}

#[cfg(feature = "cuda")]
unsafe impl<LeafData, Storage> DeviceCopy for GenericQbvh<LeafData, Storage>
where
    Storage: QbvhStorage<LeafData>,
    Storage::Nodes: DeviceCopy,
    Storage::ArrayU32: DeviceCopy,
    Storage::ArrayProxies: DeviceCopy,
{
}

#[cfg(all(feature = "std", feature = "cuda"))]
impl<LeafData: DeviceCopy> CudaQbvh<LeafData> {
    /// Returns the qbvh usable from within a CUDA kernel.
    pub fn as_device_ptr(&self) -> CudaQbvhPtr<LeafData> {
        GenericQbvh {
            root_aabb: self.root_aabb,
            nodes: self.nodes.as_device_ptr(),
            dirty_nodes: self.dirty_nodes.as_device_ptr(),
            free_list: self.free_list.as_device_ptr(),
            proxies: self.proxies.as_device_ptr(),
        }
    }
}

#[cfg(feature = "std")]
impl<LeafData: IndexedData> Qbvh<LeafData> {
    /// Initialize an empty Qbvh.
    pub fn new() -> Self {
        Qbvh {
            root_aabb: Aabb::new_invalid(),
            nodes: Vec::new(),
            dirty_nodes: Vec::new(),
            free_list: Vec::new(),
            proxies: Vec::new(),
        }
    }

    /// Converts this RAM-based qbvh to an qbvh based on CUDA memory.
    #[cfg(feature = "cuda")]
    pub fn to_cuda(&self) -> CudaResult<CudaQbvh<LeafData>>
    where
        LeafData: DeviceCopy,
    {
        Ok(CudaQbvh {
            root_aabb: self.root_aabb,
            nodes: CudaArray1::new(&self.nodes)?,
            dirty_nodes: CudaArray1::new(&self.dirty_nodes)?,
            free_list: CudaArray1::new(&self.free_list)?,
            proxies: CudaArray1::new(&self.proxies)?,
        })
    }

    /// Iterates mutably through all the leaf data in this Qbvh.
    pub fn iter_data_mut(&mut self) -> impl Iterator<Item = (NodeIndex, &mut LeafData)> {
        self.proxies.iter_mut().map(|p| (p.node, &mut p.data))
    }

    /// Iterate through all the leaf data in this Qbvh.
    pub fn iter_data(&self) -> impl Iterator<Item = (NodeIndex, &LeafData)> {
        self.proxies.iter().map(|p| (p.node, &p.data))
    }

    /// The Aabb of the given node.
    pub fn node_aabb(&self, node_id: NodeIndex) -> Option<Aabb> {
        self.nodes
            .get(node_id.index as usize)
            .map(|n| n.simd_aabb.extract(node_id.lane as usize))
    }

    /// Returns the data associated to a given leaf.
    ///
    /// Returns `None` if the provided node ID does not identify a leaf.
    pub fn leaf_data(&self, node_id: NodeIndex) -> Option<LeafData> {
        let node = self.nodes.get(node_id.index as usize)?;

        if !node.is_leaf() {
            return None;
        }

        let proxy = self
            .proxies
            .get(node.children[node_id.lane as usize] as usize)?;
        Some(proxy.data)
    }

    /// The raw nodes of this BVH.
    ///
    /// If this Qbvh isn’t empty, the first element of the returned slice is the root of the
    /// tree. The other elements are not arranged in any particular order.
    /// The more high-level traversal methods should be used instead of this.
    pub fn raw_nodes(&self) -> &[QbvhNode] {
        &self.nodes
    }

    /// The raw proxies of this BVH.
    ///
    /// If this Qbvh isn’t empty, the first element of the returned slice is the root of the
    /// tree. The other elements are not arranged in any particular order.
    /// The more high-level traversal methods should be used instead of this.
    pub fn raw_proxies(&self) -> &[QbvhProxy<LeafData>] {
        &self.proxies
    }

    /// Computes a scaled version of this Qbvh.
    ///
    /// This will apply the scale to each Aabb on this BVH.
    pub fn scaled(mut self, scale: &Vector<Real>) -> Self {
        self.root_aabb = self.root_aabb.scaled(scale);
        for node in &mut self.nodes {
            node.simd_aabb = node.simd_aabb.scaled(&Vector::splat(*scale));
        }
        self
    }
}

impl<LeafData: IndexedData, Storage: QbvhStorage<LeafData>> GenericQbvh<LeafData, Storage> {
    /// The Aabb of the root of this tree.
    pub fn root_aabb(&self) -> &Aabb {
        &self.root_aabb
    }
}

#[cfg(test)]
mod test {
    use crate::bounding_volume::Aabb;
    use crate::math::{Point, Vector};
    use crate::partitioning::Qbvh;

    #[test]
    fn multiple_identical_aabb_stack_overflow() {
        // A stack overflow was caused during the construction of the
        // Qbvh with more than four Aabb with the same center.
        let aabb = Aabb::new(Point::origin(), Vector::repeat(1.0).into());

        for k in 0u32..20 {
            let mut tree = Qbvh::new();
            tree.clear_and_rebuild((0..k).map(|i| (i, aabb)), 0.0);
        }
    }
}
