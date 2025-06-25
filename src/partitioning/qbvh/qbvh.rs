use crate::bounding_volume::{Aabb, SimdAabb};
use crate::math::{Real, Vector};
use alloc::vec::Vec;

use na::SimdValue;

#[cfg(feature = "rkyv")]
use rkyv::{bytecheck, CheckBytes};

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
        u64::MAX
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

#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Deserialize, rkyv::Serialize),
    archive(as = "Self")
)]
#[derive(Clone, Copy, Debug, Default, Eq, Hash, Ord, PartialEq, PartialOrd)]
/// The status of a QBVH node.
pub struct QbvhNodeFlags(u8);

bitflags::bitflags! {
    impl QbvhNodeFlags: u8 {
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
    /// Does the Aabb of this node needs to be updated?
    pub fn is_dirty(&self) -> bool {
        self.flags.contains(QbvhNodeFlags::DIRTY)
    }

    #[inline]
    /// Sets if the Aabb of this node needs to be updated.
    pub fn set_dirty(&mut self, dirty: bool) {
        self.flags.set(QbvhNodeFlags::DIRTY, dirty);
    }

    #[inline]
    /// Was the Aabb of this node changed since the last rebalancing?
    pub fn is_changed(&self) -> bool {
        self.flags.contains(QbvhNodeFlags::CHANGED)
    }

    #[inline]
    /// Sets if the Aabb of this node changed since the last rebalancing.
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
#[repr(C)]
#[derive(Debug, Clone)]
pub struct Qbvh<LeafData> {
    pub(super) root_aabb: Aabb,
    pub(super) nodes: Vec<QbvhNode>,
    pub(super) dirty_nodes: Vec<u32>,
    pub(super) free_list: Vec<u32>,
    pub(super) proxies: Vec<QbvhProxy<LeafData>>,
}

impl<LeafData: IndexedData> Default for Qbvh<LeafData> {
    fn default() -> Self {
        Self::new()
    }
}

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

    // TODO: support a crate like get_size2 (will require support on nalgebra too)?
    /// An approximation of the memory usage (in bytes) for this struct plus
    /// the memory it allocates dynamically.
    pub fn total_memory_size(&self) -> usize {
        size_of::<Self>() + self.heap_memory_size()
    }

    /// An approximation of the memory dynamically-allocated by this struct.
    pub fn heap_memory_size(&self) -> usize {
        let Self {
            root_aabb: _,
            nodes,
            dirty_nodes,
            free_list,
            proxies,
        } = self;
        let sz_root_aabb = 0;
        let sz_nodes = nodes.capacity() * size_of::<u32>();
        let sz_dirty_nodes = dirty_nodes.capacity() * size_of::<u32>();
        let sz_free_list = free_list.capacity() * size_of::<u32>();
        let sz_proxies = proxies.capacity() * size_of::<u32>();

        sz_root_aabb + sz_nodes + sz_dirty_nodes + sz_free_list + sz_proxies
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

impl<LeafData: IndexedData> Qbvh<LeafData> {
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
