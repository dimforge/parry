use crate::bounding_volume::{SimdAABB, AABB};
use crate::math::{Real, Vector};
use na::SimdValue;
use std::collections::VecDeque;

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

/// The index of a node part of a QBVH.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
#[cfg_attr(feature = "serde-serialize", derive(Serialize, Deserialize))]
pub struct NodeIndex {
    pub(super) index: u32, // Index of the addressed node in the `nodes` array.
    pub(super) lane: u8,   // SIMD lane of the addressed node.
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
}

/// A SIMD node of an SIMD QBVH.
///
/// This groups four nodes of the QBVH.
#[derive(Copy, Clone, Debug)]
#[cfg_attr(feature = "serde-serialize", derive(Serialize, Deserialize))]
pub struct QBVHNode {
    /// The AABBs of the qbvh nodes represented by this node.
    pub simd_aabb: SimdAABB,
    /// Index of the nodes of the 4 nodes represented by `self`.
    /// If this is a leaf, it contains the proxy ids instead.
    pub children: [u32; 4],
    /// The index of the node parent to the 4 nodes represented by `self`.
    pub parent: NodeIndex,
    /// Are the four nodes represented by `self` leaves of the `QBVH`?
    pub leaf: bool, // TODO: pack this with the NodexIndex.lane?
    pub(super) dirty: bool, // TODO: move this to a separate bitvec?
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
#[cfg_attr(feature = "serde-serialize", derive(Serialize, Deserialize))]
pub struct QBVHProxy<T> {
    pub node: NodeIndex,
    pub data: T, // The collider data. TODO: only set the collider generation here?
}

impl<T> QBVHProxy<T> {
    pub(super) fn invalid() -> Self
    where
        T: IndexedData,
    {
        Self {
            node: NodeIndex::invalid(),
            data: T::default(),
        }
    }

    pub(super) fn detached(data: T) -> Self {
        Self {
            node: NodeIndex::invalid(),
            data,
        }
    }
}

/// A quaternary bounding-volume-hierarchy.
///
/// This is a bounding-volume-hierarchy where each node has either four children or none.
#[cfg_attr(feature = "serde-serialize", derive(Serialize, Deserialize))]
#[derive(Clone, Debug)]
pub struct QBVH<T> {
    pub(super) root_aabb: AABB,
    pub(super) nodes: Vec<QBVHNode>,
    pub(super) dirty_nodes: VecDeque<u32>,
    pub(super) proxies: Vec<QBVHProxy<T>>,
}

impl<T: IndexedData> QBVH<T> {
    /// Initialize an empty QBVH.
    pub fn new() -> Self {
        QBVH {
            root_aabb: AABB::new_invalid(),
            nodes: Vec::new(),
            dirty_nodes: VecDeque::new(),
            proxies: Vec::new(),
        }
    }

    /// Iterates mutably through all the leaf data in this QBVH.
    pub fn iter_data_mut(&mut self) -> impl Iterator<Item = (NodeIndex, &mut T)> {
        self.proxies.iter_mut().map(|p| (p.node, &mut p.data))
    }

    /// Iterate through all the leaf data in this QBVH.
    pub fn iter_data(&self) -> impl Iterator<Item = (NodeIndex, &T)> {
        self.proxies.iter().map(|p| (p.node, &p.data))
    }

    /// The AABB of the root of this tree.
    pub fn root_aabb(&self) -> &AABB {
        &self.root_aabb
    }

    /// The AABB of the given node.
    pub fn node_aabb(&self, node_id: NodeIndex) -> Option<AABB> {
        self.nodes
            .get(node_id.index as usize)
            .map(|n| n.simd_aabb.extract(node_id.lane as usize))
    }

    /// Returns the data associated to a given leaf.
    ///
    /// Returns `None` if the provided node ID does not identify a leaf.
    pub fn leaf_data(&mut self, node_id: NodeIndex) -> Option<T> {
        let node = self.nodes.get(node_id.index as usize)?;

        if !node.leaf {
            return None;
        }

        let proxy = self
            .proxies
            .get(node.children[node_id.lane as usize] as usize)?;
        Some(proxy.data)
    }

    /// The raw nodes of this BVH.
    ///
    /// If this QBVH isn’t empty, the first element of the returned slice is the root of the
    /// tree. The other elements are not arranged in any particular order.
    /// The more high-level traversal methods should be used instead of this.
    pub fn raw_nodes(&self) -> &[QBVHNode] {
        &self.nodes
    }

    /// The raw proxies of this BVH.
    ///
    /// If this QBVH isn’t empty, the first element of the returned slice is the root of the
    /// tree. The other elements are not arranged in any particular order.
    /// The more high-level traversal methods should be used instead of this.
    pub fn raw_proxies(&self) -> &[QBVHProxy<T>] {
        &self.proxies
    }

    /// Computes a scaled version of this QBVH.
    ///
    /// This will apply the scale to each AABB on this BVH.
    pub fn scaled(mut self, scale: &Vector<Real>) -> Self {
        self.root_aabb = self.root_aabb.scaled(scale);
        for node in &mut self.nodes {
            node.simd_aabb = node.simd_aabb.scaled(&Vector::splat(*scale));
        }
        self
    }
}

#[cfg(test)]
mod test {
    use crate::bounding_volume::AABB;
    use crate::math::{Point, Vector};
    use crate::partitioning::QBVH;

    #[test]
    fn multiple_identical_aabb_stack_overflow() {
        // A stack overflow was caused during the construction of the
        // QBVH with more than four AABB with the same center.
        let aabb = AABB::new(Point::origin(), Vector::repeat(1.0).into());

        for k in 0u32..20 {
            let mut tree = QBVH::new();
            tree.clear_and_rebuild((0..k).map(|i| (i, aabb)), 0.0);
        }
    }
}
