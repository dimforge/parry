use super::BvhOptimizationHeapEntry;
use crate::bounding_volume::{Aabb, BoundingVolume};
use crate::math::{Point, Real, Vector};
use crate::query::{QueryOptionsNotUsed, Ray, RayCast};
use crate::utils::VecMap;
use alloc::collections::{BinaryHeap, VecDeque};
use alloc::vec::Vec;
use core::ops::{Deref, DerefMut, Index, IndexMut};

/// The strategy for one-time build of the tree.
///
/// For general-purpose usage [`BvhBuildStrategy::Binned`] is recommended. If the focus is
/// ray-casting, [`BvhBuildStrategy::Ploc`] is recommended.
#[derive(Default, Clone, Debug, Copy, PartialEq, Eq)]
pub enum BvhBuildStrategy {
    /// The tree is built using the binned strategy.
    ///
    /// This implements the strategy from "On fast Construction of SAH-based Bounding Volume Hierarchies", Ingo Ward.
    /// This is **not** parallelized yet.
    #[default]
    Binned,
    /// The tree is built using the Locally-Ordered Clustering technique.
    ///
    /// This implements the strategy from "Parallel Locally-Ordered Clustering for Bounding Volume Hierarchy Construction", Meister, Bittner.
    /// This is **not** parallelized yet.
    Ploc,
}

/// Workspace data for various operations on the tree.
///
/// This is all temporary data that can be freed at any time without affecting results.
/// The main reason to reuse the same instance of this over time is to lower costs of internal
/// allocations.
#[derive(Clone, Default)]
pub struct BvhWorkspace {
    pub(super) refit_tmp: BvhNodeVec,
    pub(super) rebuild_leaves: Vec<BvhNode>,
    pub(super) rebuild_frame_index: u32,
    pub(super) rebuild_start_index: u32,
    pub(super) optimization_roots: Vec<u32>,
    pub(super) queue: BinaryHeap<BvhOptimizationHeapEntry>,
    pub(super) dequeue: VecDeque<u32>,
    pub(super) traversal_stack: Vec<u32>,
}

/// A piece of data packing state flags as well as leaf counts for a BVH tree node.
#[derive(Default, Copy, Clone, Debug)]
#[cfg_attr(feature = "serde-serialize", derive(Serialize, Deserialize))]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Deserialize, rkyv::Serialize),
    archive(check_bytes)
)]
#[repr(transparent)]
pub struct BvhNodeData(u32);
const CHANGED: u32 = 0b01;
const CHANGE_PENDING: u32 = 0b11;

impl BvhNodeData {
    #[inline(always)]
    pub(super) fn with_leaf_count_and_pending_change(leaf_count: u32) -> Self {
        Self(leaf_count | (CHANGE_PENDING << 30))
    }

    #[inline(always)]
    pub(super) fn leaf_count(self) -> u32 {
        self.0 & 0x3fff_ffff
    }

    #[inline(always)]
    pub(super) fn is_changed(self) -> bool {
        self.0 >> 30 == CHANGED
    }

    #[inline(always)]
    pub(super) fn is_change_pending(self) -> bool {
        self.0 >> 30 == CHANGE_PENDING
    }

    #[inline(always)]
    pub(super) fn add_leaf_count(&mut self, added: u32) {
        self.0 += added;
    }

    #[inline(always)]
    pub(super) fn set_change_pending(&mut self) {
        self.0 |= CHANGE_PENDING << 30;
    }

    #[inline(always)]
    pub(super) fn resolve_pending_change(&mut self) {
        if self.is_change_pending() {
            *self = Self((self.0 & 0x3fff_ffff) | (CHANGED << 30));
        } else {
            *self = Self(self.0 & 0x3fff_ffff);
        }
    }

    pub(super) fn merged(self, other: Self) -> Self {
        let leaf_count = self.leaf_count() + other.leaf_count();
        let changed = (self.0 >> 30) | (other.0 >> 30);
        Self(leaf_count | changed << 30)
    }
}

/// A pair of tree nodes.
///
/// Both `left` and `right` are guaranteed to be valid except for the only special-case where the
/// tree contains only a single leaf, in which case only `left` is valid. But in every other
/// cases where the tree contains at least 2 leaves, booth `left` and `right` are guaranteed
/// to be valid.
#[derive(Copy, Clone, Debug)]
#[cfg_attr(feature = "serde-serialize", derive(Serialize, Deserialize))]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Deserialize, rkyv::Serialize),
    archive(check_bytes)
)]
#[repr(C)] // SAFETY: needed to ensure SIMD aabb checks rely on the layout.
pub struct BvhNodeWide {
    pub(super) left: BvhNode,
    pub(super) right: BvhNode,
}

impl BvhNodeWide {
    #[inline(always)]
    pub fn zeros() -> Self {
        Self {
            left: BvhNode::zeros(),
            right: BvhNode::zeros(),
        }
    }

    /// The two nodes in `self` seen as an array.
    ///
    /// Useful for accessing the nodes by index.
    #[inline(always)]
    pub fn as_array(&self) -> [&BvhNode; 2] {
        [&self.left, &self.right]
    }

    /// The two nodes in `self` seen as an array of mutable references.
    ///
    /// Useful for accessing the nodes mutable by index.
    #[inline(always)]
    pub fn as_array_mut(&mut self) -> [&mut BvhNode; 2] {
        [&mut self.left, &mut self.right]
    }

    /// Merges both nodes contained by `self` to form its parent.
    pub fn merged(&self, my_id: u32) -> BvhNode {
        self.left.merged(&self.right, my_id)
    }

    /// The sum of leaves contained by both nodes in `self`.
    pub fn leaf_count(&self) -> u32 {
        self.left.leaf_count() + self.right.leaf_count()
    }
}

#[repr(C)] // SAFETY: needed to ensure SIMD aabb checks rely on the layout.
#[cfg(all(feature = "simd-is-enabled", feature = "dim3", feature = "f32"))]
pub(super) struct BvhNodeSimd {
    mins: glam::Vec3A,
    maxs: glam::Vec3A,
}

// SAFETY: compile-time assertions to ensure we can transmute between `BvhNode` and `BvhNodeSimd`.
#[cfg(all(feature = "simd-is-enabled", feature = "dim3", feature = "f32"))]
static_assertions::assert_eq_align!(BvhNode, BvhNodeSimd);
#[cfg(all(feature = "simd-is-enabled", feature = "dim3", feature = "f32"))]
static_assertions::assert_eq_size!(BvhNode, BvhNodeSimd);

/// The note (internal or leaf) of a BVH.
#[derive(Copy, Clone, Debug)]
#[repr(C)] // SAFETY: needed to ensure SIMD aabb checks rely on the layout.
#[cfg_attr(all(feature = "f32", feature = "dim3"), repr(align(16)))]
#[cfg_attr(feature = "serde-serialize", derive(Serialize, Deserialize))]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Deserialize, rkyv::Serialize),
    archive(check_bytes)
)]
pub struct BvhNode {
    /// Mins coordinates of the node’s bounding volume.
    pub(super) mins: Point<Real>,
    /// Children of this node. A node has either 0 (i.e. it’s a leaf) or 2 children.
    ///
    /// If [`Self::leaf_count`] is 1, then the node has 0 children and is a leaf.
    pub(super) children: u32,
    /// Maxs coordinates of this node’s bounding volume.
    pub(super) maxs: Point<Real>,
    /// Packed data associated to this node (leaf count and flags).
    pub(super) data: BvhNodeData,
}

impl BvhNode {
    #[inline(always)]
    pub(super) fn zeros() -> Self {
        Self {
            mins: Point::origin(),
            children: 0,
            maxs: Point::origin(),
            data: BvhNodeData(0),
        }
    }

    /// Initialized a leaf.
    #[inline(always)]
    pub fn leaf(aabb: Aabb, leaf_data: u32) -> BvhNode {
        Self {
            mins: aabb.mins,
            maxs: aabb.maxs,
            children: leaf_data,
            data: BvhNodeData::with_leaf_count_and_pending_change(1),
        }
    }

    /// If this node is a leaf, returns its associated index provided at construction time.
    #[inline(always)]
    pub fn leaf_data(&self) -> Option<u32> {
        self.is_leaf().then_some(self.children)
    }

    /// Is this node a leaf?
    #[inline(always)]
    pub fn is_leaf(&self) -> bool {
        self.leaf_count() == 1
    }

    #[inline(always)]
    pub(super) fn leaf_count(&self) -> u32 {
        self.data.leaf_count()
    }

    #[inline(always)]
    #[cfg(all(feature = "simd-is-enabled", feature = "dim3", feature = "f32"))]
    pub(super) fn as_simd(&self) -> &BvhNodeSimd {
        // SAFETY: BvhNode is declared with the alignment
        //         and size of two SimdReal.
        unsafe { core::mem::transmute(self) }
    }

    #[inline(always)]
    pub(super) fn merged(&self, other: &Self, children: u32) -> Self {
        // TODO PERF: simd optimizations?
        Self {
            mins: self.mins.inf(&other.mins),
            children,
            maxs: self.maxs.sup(&other.maxs),
            data: self.data.merged(other.data),
        }
    }

    /// The min corner of this node’s AABB.
    #[inline]
    pub fn mins(&self) -> Point<Real> {
        self.mins
    }

    /// The max corner of this node’s AABB.
    #[inline]
    pub fn maxs(&self) -> Point<Real> {
        self.maxs
    }

    /// This node’s AABB.
    #[inline]
    pub fn aabb(&self) -> Aabb {
        Aabb {
            mins: self.mins,
            maxs: self.maxs,
        }
    }

    /// The center of this node’s AABB.
    #[inline]
    pub fn center(&self) -> Point<Real> {
        na::center(&self.mins, &self.maxs)
    }

    /// Return `true` if this node has been marked as changed.
    #[inline(always)]
    pub fn is_changed(&self) -> bool {
        self.data.is_changed()
    }

    /// Applies a scale factor to this node’s AABB.
    #[inline]
    pub fn scale(&mut self, scale: &Vector<Real>) {
        self.mins.coords.component_mul_assign(scale);
        self.maxs.coords.component_mul_assign(scale);
    }

    /// Calculates the volume of the AABB of `self`.
    #[inline]
    pub fn volume(&self) -> Real {
        // TODO PERF: simd optimizations?
        let extents = self.maxs - self.mins;
        #[cfg(feature = "dim2")]
        return extents.x * extents.y;
        #[cfg(feature = "dim3")]
        return extents.x * extents.y * extents.z;
    }

    /// Calculates the volume of the AABB resulting from the merge of `self` and `other`.
    pub fn merged_volume(&self, other: &Self) -> Real {
        // TODO PERF: simd optimizations?
        let mins = self.mins.inf(&other.mins);
        let maxs = self.maxs.sup(&other.maxs);
        let extents = maxs - mins;

        #[cfg(feature = "dim2")]
        return extents.x * extents.y;
        #[cfg(feature = "dim3")]
        return extents.x * extents.y * extents.z;
    }

    /// Checks if the AABB of `self` intersects the `other` node’s AABB.
    #[cfg(not(all(feature = "simd-is-enabled", feature = "dim3", feature = "f32")))]
    pub fn intersects(&self, other: &Self) -> bool {
        na::partial_le(&self.mins, &other.maxs) && na::partial_ge(&self.maxs, &other.mins)
    }

    /// Checks if the AABB of `self` intersects the `other` node’s AABB.
    #[cfg(all(feature = "simd-is-enabled", feature = "dim3", feature = "f32"))]
    pub fn intersects(&self, other: &Self) -> bool {
        let simd_self = self.as_simd();
        let simd_other = other.as_simd();
        (simd_self.mins.cmple(simd_other.maxs) & simd_self.maxs.cmpge(simd_other.mins)).all()
    }

    /// Checks if the AABB of `self` fully encloses the `other` node’s AABB.
    #[cfg(not(all(feature = "simd-is-enabled", feature = "dim3", feature = "f32")))]
    pub fn contains(&self, other: &Self) -> bool {
        na::partial_le(&self.mins, &other.mins) && na::partial_ge(&self.maxs, &other.maxs)
    }

    /// Checks if the AABB of `self` fully encloses the `other` node’s AABB.
    #[cfg(all(feature = "simd-is-enabled", feature = "dim3", feature = "f32"))]
    pub fn contains(&self, other: &Self) -> bool {
        let simd_self = self.as_simd();
        let simd_other = other.as_simd();
        (simd_self.mins.cmple(simd_other.mins) & simd_self.maxs.cmpge(simd_other.maxs)).all()
    }

    /// Checks if the AABB of `self` fully encloses the `other` AABB.
    pub fn contains_aabb(&self, other: &Aabb) -> bool {
        // TODO PERF: simd optimizations?
        na::partial_le(&self.mins, &other.mins) && na::partial_ge(&self.maxs, &other.maxs)
    }

    /// Casts a ray on this AABB.
    ///
    /// Returns `Real::MAX` if there is no hit.
    pub fn cast_ray(&self, ray: &Ray, max_toi: Real) -> Real {
        self.aabb()
            .cast_local_ray(ray, max_toi, true, &QueryOptionsNotUsed)
            .unwrap_or(Real::MAX)
    }

    /// Casts a ray on this AABB, with SIMD optimizations.
    ///
    /// Returns `Real::MAX` if there is no hit.
    #[cfg(all(feature = "simd-is-enabled", feature = "dim3", feature = "f32"))]
    pub(super) fn cast_inv_ray_simd(&self, ray: &super::bvh_queries::SimdInvRay) -> f32 {
        let simd_self = self.as_simd();
        let t1 = (simd_self.mins - ray.origin) * ray.inv_dir;
        let t2 = (simd_self.maxs - ray.origin) * ray.inv_dir;

        let tmin = t1.min(t2);
        let tmax = t1.max(t2);
        // let tmin = tmin.as_array_ref();
        // let tmax = tmax.as_array_ref();
        let tmin_n = tmin.max_element(); // tmin[0].max(tmin[1].max(tmin[2]));
        let tmax_n = tmax.min_element(); // tmax[0].min(tmax[1].min(tmax[2]));

        if tmax_n >= tmin_n && tmax_n >= 0.0 {
            tmin_n
        } else {
            f32::MAX
        }
    }
}

/// An index identifying a single BVH tree node.
#[derive(Copy, Clone, Debug, Default, PartialEq, Eq)]
#[cfg_attr(feature = "serde-serialize", derive(Serialize, Deserialize))]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Deserialize, rkyv::Serialize),
    archive(check_bytes)
)]
pub struct BvhNodeIndex(pub usize);

impl BvhNodeIndex {
    pub(super) const LEFT: bool = false;
    pub(super) const RIGHT: bool = true;

    /// Decomposes the index into the `BvhNodeWide` and the boolean indicating if
    /// it is `BvhNodeWide::left` (false) or `BvhNodeWide::right` (true).
    #[inline]
    pub fn decompose(self) -> (usize, bool) {
        (self.0 >> 1, (self.0 & 0b01) != 0)
    }

    /// The BVH index for the sibling of `self`.
    ///
    /// For example if `self` identifies the left child of a node, then its sibling is the right
    /// child of the same node.
    #[inline]
    pub fn sibling(self) -> Self {
        // Just flip the first bit to switch between left and right child.
        Self(self.0 ^ 0b01)
    }

    /// The BVH node index stored as the left child of the `id`-th 2-wide node of this tree.
    #[inline]
    pub fn left(id: u32) -> Self {
        Self::new(id, Self::LEFT)
    }

    /// The BVH node index stored as the right child of the `id`-th 2-wide node of this tree.
    #[inline]
    pub fn right(id: u32) -> Self {
        Self::new(id, Self::RIGHT)
    }

    /// The BVH node index stored as the left (`is_right = false`) or right (`is_right = true`)
    /// child of the `id`-th 2-wide node of this tree.
    #[inline]
    pub fn new(id: u32, is_right: bool) -> Self {
        Self(((id as usize) << 1) | (is_right as usize))
    }
}

#[derive(Clone, Debug, Default)]
#[cfg_attr(feature = "serde-serialize", derive(Serialize, Deserialize))]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Deserialize, rkyv::Serialize),
    archive(check_bytes)
)]
pub(crate) struct BvhNodeVec(pub(crate) Vec<BvhNodeWide>);

impl Deref for BvhNodeVec {
    type Target = Vec<BvhNodeWide>;
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl DerefMut for BvhNodeVec {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl Index<usize> for BvhNodeVec {
    type Output = BvhNodeWide;

    #[inline(always)]
    fn index(&self, index: usize) -> &Self::Output {
        &self.0[index]
    }
}

impl IndexMut<usize> for BvhNodeVec {
    #[inline(always)]
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.0[index]
    }
}

impl Index<BvhNodeIndex> for BvhNodeVec {
    type Output = BvhNode;

    #[inline(always)]
    fn index(&self, index: BvhNodeIndex) -> &Self::Output {
        self.0[index.0 >> 1].as_array()[index.0 & 1]
    }
}

impl IndexMut<BvhNodeIndex> for BvhNodeVec {
    #[inline(always)]
    fn index_mut(&mut self, index: BvhNodeIndex) -> &mut Self::Output {
        self.0[index.0 >> 1].as_array_mut()[index.0 & 1]
    }
}

/// A Bounding Volume Hierarchy designed for spatial queries and physics engine broad-phases.
#[derive(Clone, Debug, Default)]
#[cfg_attr(feature = "serde-serialize", derive(Serialize, Deserialize))]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Deserialize, rkyv::Serialize),
    archive(check_bytes)
)]
pub struct Bvh {
    pub(super) nodes: BvhNodeVec,
    // Parent indices for elements in `nodes`.
    // We don’t store this in `Self::nodes` since it’s only useful for node removal.
    pub(super) parents: Vec<BvhNodeIndex>,
    pub(super) leaf_node_indices: VecMap<BvhNodeIndex>,
}

impl Bvh {
    /// An empty BVH.
    pub fn new() -> Self {
        Self::default()
    }

    /// Creates a new BVH with a slice of AABBs.
    ///
    /// Each leaf will be associated an index equal to its position into the slice. For example,
    /// the AABB `leaves[42]` is associated to the leaf with index 42.
    pub fn from_leaves(strategy: BvhBuildStrategy, leaves: &[Aabb]) -> Self {
        Self::from_iter(strategy, leaves.iter().copied().enumerate())
    }

    /// Creates a new BVH with leaves given by an iterator.
    ///
    /// The iterator yields leaf index and aabbs. The leaf indices will then be read back
    /// by various methods like tree traversals or leaf iterations.
    ///
    /// Note that the indices are stored internally as `u32`. The iterator expects `usize`
    /// for convenience (so that iterators built with `.enumerate()` can be used directly
    /// without an additional cast of the `usize` index to `u32`).
    pub fn from_iter<It>(strategy: BvhBuildStrategy, leaves: It) -> Self
    where
        It: IntoIterator<Item = (usize, Aabb)>,
    {
        let leaves = leaves.into_iter();
        let (capacity_lo, capacity_up) = leaves.size_hint();
        let capacity = capacity_up.unwrap_or(capacity_lo);

        let mut result = Self::new();
        let mut workspace = BvhWorkspace::default();
        workspace.rebuild_leaves.reserve(capacity);
        result.leaf_node_indices.reserve_len(capacity);

        for (leaf_id, leaf_aabb) in leaves {
            workspace
                .rebuild_leaves
                .push(BvhNode::leaf(leaf_aabb, leaf_id as u32));
            let _ = result
                .leaf_node_indices
                .insert(leaf_id, BvhNodeIndex::default());
        }

        // Handle special cases that don’t play well with the rebuilds.
        match workspace.rebuild_leaves.len() {
            0 => {}
            1 => {
                result.nodes.push(BvhNodeWide {
                    left: workspace.rebuild_leaves[0],
                    right: BvhNode::zeros(),
                });
                result.parents.push(BvhNodeIndex::default());
            }
            2 => {
                result.nodes.push(BvhNodeWide {
                    left: workspace.rebuild_leaves[0],
                    right: workspace.rebuild_leaves[1],
                });
                result.parents.push(BvhNodeIndex::default());
            }
            _ => {
                result.nodes.reserve(capacity);
                result.parents.reserve(capacity);
                result.parents.clear();
                result.nodes.push(BvhNodeWide::zeros());
                result.parents.push(BvhNodeIndex::default());

                match strategy {
                    BvhBuildStrategy::Ploc => {
                        result.rebuild_range_ploc(0, &mut workspace.rebuild_leaves)
                    }
                    BvhBuildStrategy::Binned => {
                        result.rebuild_range_binned(0, &mut workspace.rebuild_leaves)
                    }
                }

                // Layout in depth-first order.
                result.refit(&mut workspace);
            }
        }

        result
    }

    /// The AABB bounding everything contained by this BVH.
    pub fn root_aabb(&self) -> Aabb {
        match self.leaf_count() {
            0 => Aabb::new_invalid(),
            1 => self.nodes[0].left.aabb(),
            _ => self.nodes[0]
                .left
                .aabb()
                .merged(&self.nodes[0].right.aabb()),
        }
    }

    /// Scales this tree’s geometry by the given factor.
    ///
    /// This is only valid if all elements of `scale` are positive.
    pub fn scale(&mut self, scale: &Vector<Real>) {
        for node in self.nodes.0.iter_mut() {
            node.left.scale(scale);
            node.right.scale(scale);
        }
    }

    /// Does this tree not contain any leaf?
    pub fn is_empty(&self) -> bool {
        self.nodes.is_empty()
    }

    ///  Reference to the leaf associated to the given index at construction time.
    pub fn leaf_node(&self, leaf_key: u32) -> Option<&BvhNode> {
        let idx = self.leaf_node_indices.get(leaf_key as usize)?;
        Some(&self.nodes[*idx])
    }

    /// An approximation of the memory usage (in bytes) for this struct plus
    /// the memory it allocates dynamically.
    pub fn total_memory_size(&self) -> usize {
        size_of::<Self>() + self.heap_memory_size()
    }

    /// An approximation of the memory dynamically-allocated by this struct.
    pub fn heap_memory_size(&self) -> usize {
        let Self {
            nodes,
            parents,
            leaf_node_indices,
        } = self;
        nodes.capacity() * size_of::<BvhNodeWide>()
            + parents.capacity() * size_of::<BvhNodeIndex>()
            + leaf_node_indices.capacity() * size_of::<BvhNodeIndex>()
    }

    /// The depth of the sub-tree rooted at the node with index `node_id`.
    ///
    /// Set `node_id` to 0 to get the depth of the whole tree.
    pub fn subtree_depth(&self, node_id: u32) -> u32 {
        if node_id == 0 && self.nodes.is_empty() {
            return 0;
        } else if node_id == 0 && self.nodes.len() == 1 {
            return 1 + (self.nodes[0].right.leaf_count() != 0) as u32;
        }

        let node = &self.nodes[node_id as usize];

        let left_depth = if node.left.is_leaf() {
            1
        } else {
            self.subtree_depth(node.left.children)
        };

        let right_depth = if node.right.is_leaf() {
            1
        } else {
            self.subtree_depth(node.right.children)
        };

        left_depth.max(right_depth) + 1
    }

    /// The number of leaves of this tree.
    pub fn leaf_count(&self) -> u32 {
        if self.nodes.is_empty() {
            0
        } else {
            self.nodes[0].leaf_count()
        }
    }

    /// Deletes the leaf at the given index.
    ///
    /// The operation is `O(h)` where `h` is the tree height (which, if
    /// optimized should be close to `n log(n)` where `n` is the leaf count). It will update the
    /// leaf counts and AABBs of all ancestors of the removed leaf.
    // TODO: should we make a version that doesn’t traverse the parents?
    //       If we do, we must be very careful that the leaf counts that become
    //       invalid don’t break other algorithm… (and, in particular, the root
    //       special case that checks if its right element has 0 leaf count).
    pub fn remove(&mut self, leaf_index: u32) {
        if let Some(node_index) = self.leaf_node_indices.remove(leaf_index as usize) {
            if self.leaf_node_indices.is_empty() {
                // We deleted the last leaf! Remove the root.
                self.nodes.clear();
                self.parents.clear();
                return;
            }

            let sibling = node_index.sibling();
            let (wide_node_index, is_right) = node_index.decompose();

            if wide_node_index == 0 {
                if self.nodes[sibling].is_leaf() {
                    // If the sibling is a leaf, we end up with a partial root.
                    // There is no parent pointer to update.
                    if !is_right {
                        // We remove the left leaf. Move the right leaf in its place.
                        let moved_index = self.nodes[0].right.children;
                        self.nodes[0].left = self.nodes[0].right;
                        self.leaf_node_indices[moved_index as usize] = BvhNodeIndex::left(0);
                    }

                    // Now we can just clear the right leaf.
                    self.nodes[0].right = BvhNode::zeros();
                } else {
                    // The sibling isn’t a leaf. It becomes the new root at index 0.
                    self.nodes[0] = self.nodes[sibling.decompose().0];
                    // Both parent pointers need to be updated since both nodes moved to the root.
                    let new_root = &mut self.nodes[0];
                    if new_root.left.is_leaf() {
                        self.leaf_node_indices[new_root.left.children as usize] =
                            BvhNodeIndex::left(0);
                    } else {
                        self.parents[new_root.left.children as usize] = BvhNodeIndex::left(0);
                    }
                    if new_root.right.is_leaf() {
                        self.leaf_node_indices[new_root.right.children as usize] =
                            BvhNodeIndex::right(0);
                    } else {
                        self.parents[new_root.right.children as usize] = BvhNodeIndex::right(0);
                    }
                }
            } else {
                // The sibling moves to the parent. The affected wide node is no longer accessible,
                // but we can just leave it there, it will get cleaned up at the next refit.
                let parent = self.parents[wide_node_index];
                let sibling = &self.nodes[sibling];

                if sibling.is_leaf() {
                    self.leaf_node_indices[sibling.children as usize] = parent;
                } else {
                    self.parents[sibling.children as usize] = parent;
                }

                self.nodes[parent] = *sibling;

                // TODO: we could use that propagation as an opportunity to
                //       apply some rotations?
                let mut curr = parent.decompose().0;
                while curr != 0 {
                    let parent = self.parents[curr];
                    self.nodes[parent] = self.nodes[curr].merged(curr as u32);
                    curr = parent.decompose().0;
                }
            }
        }
    }

    // pub fn quality_metric(&self) -> Real {
    //     let mut metric = 0.0;
    //     for i in 0..self.nodes.len() {
    //         if !self.nodes[i].is_leaf() {
    //             metric += self.sah_cost(i);
    //         }
    //     }
    //     metric
    // }
}
