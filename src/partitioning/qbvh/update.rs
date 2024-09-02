use crate::bounding_volume::{Aabb, BoundingVolume, SimdAabb};
#[cfg(feature = "dim3")]
use crate::math::Vector;
use crate::math::{Point, Real};
use crate::partitioning::{CenterDataSplitter, QbvhProxy};
use crate::simd::{SimdReal, SIMD_WIDTH};
use simba::simd::{SimdBool, SimdValue};

use super::{IndexedData, NodeIndex, Qbvh, QbvhNode, QbvhNodeFlags};

#[cfg(test)]
mod tests;

#[derive(Clone, Default)]
/// Workspace for QBVH update.
///
/// Re-using the same workspace for multiple QBVH modification isn’t mandatory, but it improves
/// performances by avoiding useless allocation.
pub struct QbvhUpdateWorkspace {
    stack: Vec<(u32, u8)>,
    dirty_parent_nodes: Vec<u32>,
    // For rebalancing.
    to_sort: Vec<usize>,
    orig_ids: Vec<u32>,
    is_leaf: Vec<bool>,
    aabbs: Vec<Aabb>,
}

impl QbvhUpdateWorkspace {
    fn clear(&mut self) {
        self.stack.clear();
        self.dirty_parent_nodes.clear();
        self.to_sort.clear();
        self.orig_ids.clear();
        self.is_leaf.clear();
        self.aabbs.clear();
    }
}

impl<LeafData: IndexedData> Qbvh<LeafData> {
    /// Immediately remove a leaf from this QBVH.
    pub fn remove(&mut self, data: LeafData) -> Option<LeafData> {
        let id = data.index();
        let proxy = self.proxies.get_mut(id)?;
        let node = self.nodes.get_mut(proxy.node.index as usize)?;
        node.children[proxy.node.lane as usize] = u32::MAX;

        if !node.is_dirty() {
            node.set_dirty(true);
            self.dirty_nodes.push(proxy.node.index);
        }

        *proxy = QbvhProxy::invalid();

        Some(data)
    }

    /// Prepare a new leaf for insertion into this QBVH (or for update if it already exists).
    ///
    /// The insertion or update will be completely valid only after the next call to [`Qbvh::refit`].
    pub fn pre_update_or_insert(&mut self, data: LeafData) {
        if self.nodes.is_empty() {
            let mut root = QbvhNode::empty();
            root.children = [1, u32::MAX, u32::MAX, u32::MAX];
            self.nodes
                .extend_from_slice(&[root, QbvhNode::empty_leaf_with_parent(NodeIndex::new(0, 0))]);
        }

        let proxy_id = data.index();

        if self.proxies.len() <= proxy_id {
            self.proxies.resize(proxy_id + 1, QbvhProxy::invalid());
        }

        let proxy = &mut self.proxies[proxy_id];
        proxy.data = data;

        if proxy.is_detached() {
            // Attach the proxy into one of the leaf attached to the root, if there is room.
            for ii in 0..SIMD_WIDTH {
                let mut child = self.nodes[0].children[ii];

                if child == u32::MAX {
                    // Missing node, create it.
                    child = self.nodes.len() as u32;
                    self.nodes
                        .push(QbvhNode::empty_leaf_with_parent(NodeIndex::new(
                            0, ii as u8,
                        )));
                    self.nodes[0].children[ii] = child;
                }

                let child_node = &mut self.nodes[child as usize];
                // skip non-leaf node that is attached to the root
                if !child_node.is_leaf() {
                    continue;
                }

                // Insert into this node if there is room.
                for kk in 0..SIMD_WIDTH {
                    if child_node.children[kk] == u32::MAX {
                        child_node.children[kk] = proxy_id as u32;
                        proxy.node = NodeIndex::new(child, kk as u8);

                        if !child_node.is_dirty() {
                            self.dirty_nodes.push(child);
                            child_node.set_dirty(true);
                        }

                        return;
                    }
                }
            }

            // If we reached this point, we didn’t find room for our
            // new proxy. Create a new root to make more room.
            let mut old_root = self.nodes[0];
            old_root.parent = NodeIndex::new(0, 0);

            for child_id in old_root.children {
                let parent_index = self.nodes.len() as u32;

                if let Some(child_node) = self.nodes.get_mut(child_id as usize) {
                    child_node.parent.index = parent_index;
                }
            }

            let new_leaf_node_id = self.nodes.len() as u32 + 1;
            let mut new_leaf_node = QbvhNode::empty_leaf_with_parent(NodeIndex::new(0, 1));
            new_leaf_node.children[0] = proxy_id as u32;
            // The first new node contains the (unmodified) old root, the second
            // new node contains the newly added proxy.
            self.dirty_nodes.push(new_leaf_node_id);
            new_leaf_node.set_dirty(true);
            proxy.node = NodeIndex::new(new_leaf_node_id, 0);

            let new_root_children = [
                self.nodes.len() as u32,
                new_leaf_node_id,
                u32::MAX,
                u32::MAX,
            ];

            self.nodes.extend_from_slice(&[old_root, new_leaf_node]);
            self.nodes[0].children = new_root_children;
        } else {
            let node = &mut self.nodes[proxy.node.index as usize];
            if !node.is_dirty() {
                node.set_dirty(true);
                self.dirty_nodes.push(proxy.node.index);
            }
        }
    }

    #[allow(dead_code)] // Not sure yet if we want to keep this.
    pub(crate) fn clear_changed_flag(&mut self, workspace: &mut QbvhUpdateWorkspace) {
        if self.nodes.is_empty() {
            return;
        }

        workspace.clear();
        workspace.stack.push((0, 0));

        while let Some((id, _)) = workspace.stack.pop() {
            let node = &mut self.nodes[id as usize];

            if node.is_changed() {
                node.set_changed(false);

                for child in node.children {
                    if child < self.nodes.len() as u32 {
                        workspace.stack.push((child, 0));
                    }
                }
            }
        }
    }

    #[doc(hidden)]
    pub fn check_topology(&self, check_aabbs: bool, aabb_builder: impl Fn(&LeafData) -> Aabb) {
        if self.nodes.is_empty() {
            return;
        }

        let mut stack = vec![];
        stack.push(0);

        let mut node_id_found = vec![false; self.nodes.len()];
        let mut proxy_id_found = vec![false; self.proxies.len()];

        while let Some(id) = stack.pop() {
            let node = &self.nodes[id as usize];

            // No duplicate references to internal nodes.
            assert!(!node_id_found[id as usize]);
            node_id_found[id as usize] = true;

            if id != 0 {
                let parent = &self.nodes[node.parent.index as usize];

                // Check parent <-> children indices.
                assert!(!parent.is_leaf());
                assert_eq!(parent.children[node.parent.lane as usize], id);

                if check_aabbs {
                    // Make sure the parent AABB contains its child AABB.
                    assert!(
                        parent
                            .simd_aabb
                            .extract(node.parent.lane as usize)
                            .contains(&node.simd_aabb.to_merged_aabb()),
                        "failed for {id}"
                    );
                }

                // If this node is changed, its parent is changed too.
                if node.is_changed() {
                    assert!(parent.is_changed());
                }
            }

            if node.is_leaf() {
                for ii in 0..SIMD_WIDTH {
                    let proxy_id = node.children[ii];
                    if proxy_id != u32::MAX {
                        // No duplicate references to proxies.
                        assert!(!proxy_id_found[proxy_id as usize]);
                        proxy_id_found[proxy_id as usize] = true;

                        if check_aabbs {
                            // Proxy AABB is correct.
                            let aabb = node.simd_aabb.extract(ii);
                            assert!(
                                aabb.contains(&aabb_builder(&self.proxies[proxy_id as usize].data))
                            );
                        }

                        // Proxy node reference is correct.
                        assert_eq!(
                            self.proxies[proxy_id as usize].node,
                            NodeIndex::new(id, ii as u8),
                            "Failed {proxy_id}"
                        );
                    }
                }
            } else {
                for child in node.children {
                    if child != u32::MAX {
                        stack.push(child);
                    }
                }
            }
        }

        // Each proxy was visited once.
        assert_eq!(
            proxy_id_found.iter().filter(|found| **found).count(),
            self.proxies.iter().filter(|p| !p.is_detached()).count()
        );
    }

    /// Update all the nodes that have been marked as dirty by [`Qbvh::pre_update_or_insert`],
    /// and [`Qbvh::remove`].
    ///
    /// This will not alter the topology of this `Qbvh`.
    pub fn refit<F>(
        &mut self,
        margin: Real,
        workspace: &mut QbvhUpdateWorkspace,
        aabb_builder: F,
    ) -> usize
    where
        F: Fn(&LeafData) -> Aabb,
    {
        // Loop on the dirty leaves.
        workspace.clear();
        let margin = SimdReal::splat(margin);
        let mut first_iter = true;
        let mut num_changed = 0;

        while !self.dirty_nodes.is_empty() {
            while let Some(id) = self.dirty_nodes.pop() {
                // NOTE: this will cover the case where we reach the root of the tree.
                if let Some(node) = self.nodes.get(id as usize) {
                    // Compute the new aabb.
                    let mut new_aabbs = [Aabb::new_invalid(); SIMD_WIDTH];
                    for (child_id, new_aabb) in node.children.iter().zip(new_aabbs.iter_mut()) {
                        if node.is_leaf() {
                            // We are in a leaf: compute the Aabbs.
                            if let Some(proxy) = self.proxies.get(*child_id as usize) {
                                *new_aabb = aabb_builder(&proxy.data);
                            }
                        } else {
                            // We are in an internal node: compute the children's Aabbs.
                            if let Some(node) = self.nodes.get(*child_id as usize) {
                                *new_aabb = node.simd_aabb.to_merged_aabb();
                            }
                        }
                    }

                    let node = &mut self.nodes[id as usize];
                    let new_simd_aabb = SimdAabb::from(new_aabbs);
                    node.set_dirty(false);

                    if !first_iter || !node.simd_aabb.contains(&new_simd_aabb).all() {
                        node.set_changed(true);
                        node.simd_aabb = new_simd_aabb;
                        node.simd_aabb.loosen(margin);
                        let parent_id = node.parent.index;
                        num_changed += 1;

                        if let Some(parent) = self.nodes.get_mut(parent_id as usize) {
                            if !parent.is_dirty() {
                                workspace.dirty_parent_nodes.push(parent_id);
                                parent.set_dirty(true);
                            }
                        }
                    }
                }
            }

            first_iter = false;
            std::mem::swap(&mut self.dirty_nodes, &mut workspace.dirty_parent_nodes);
        }

        num_changed
    }

    /// Rebalances the `Qbvh` tree.
    ///
    /// This will modify the topology of this tree. This assumes that the leaf AABBs have
    /// already been updated with [`Qbvh::refit`].
    pub fn rebalance(&mut self, margin: Real, workspace: &mut QbvhUpdateWorkspace) {
        if self.nodes.is_empty() {
            return;
        }

        workspace.clear();

        const MIN_CHANGED_DEPTH: u8 = 5; // TODO: find a good value

        // PERF: if we have modifications past this depth, the QBVH has become very
        //       unbalanced and a full rebuild is warranted. Note that for a perfectly
        //       balanced tree to reach this threshold, it would have to contain
        //       at least 4^15 = 1.073.741.824 leaves (i.e. very unlikely in most practical
        //       use-cases).
        // TODO: work on a way to reduce the risks of imbalance.
        const FULL_REBUILD_DEPTH: u8 = 15;

        for root_child in self.nodes[0].children {
            if root_child < self.nodes.len() as u32 {
                workspace.stack.push((root_child, 1));
            }
        }

        // Collect empty slots and pseudo-leaves to sort.
        let mut force_full_rebuild = false;
        while let Some((id, depth)) = workspace.stack.pop() {
            if depth > FULL_REBUILD_DEPTH {
                force_full_rebuild = true;
                break;
            }

            let node = &mut self.nodes[id as usize];

            if node.is_leaf() {
                self.free_list.push(id);
                for ii in 0..SIMD_WIDTH {
                    let proxy_id = node.children[ii];
                    if proxy_id < self.proxies.len() as u32 {
                        workspace.orig_ids.push(proxy_id);
                        workspace.aabbs.push(node.simd_aabb.extract(ii));
                        workspace.is_leaf.push(true);
                    }
                }
            } else if node.is_changed() || depth < MIN_CHANGED_DEPTH {
                self.free_list.push(id);

                for child in node.children {
                    if child < self.nodes.len() as u32 {
                        workspace.stack.push((child, depth + 1));
                    }
                }
            } else {
                workspace.orig_ids.push(id);
                workspace.aabbs.push(node.simd_aabb.to_merged_aabb());
                workspace.is_leaf.push(false);
            }
        }

        if force_full_rebuild {
            workspace.clear();
            let all_leaves: Vec<_> = self
                .proxies
                .iter()
                .filter_map(|proxy| self.node_aabb(proxy.node).map(|aabb| (proxy.data, aabb)))
                .collect();

            self.clear_and_rebuild(all_leaves.into_iter(), 0.0);
            return;
        }

        workspace.to_sort.extend(0..workspace.orig_ids.len());
        let root_id = NodeIndex::new(0, 0);

        let mut indices = std::mem::take(&mut workspace.to_sort);
        let (id, aabb) = self.do_recurse_rebalance(&mut indices, workspace, root_id, margin);
        workspace.to_sort = indices;

        self.root_aabb = aabb;

        self.nodes[0] = QbvhNode {
            simd_aabb: SimdAabb::from([
                aabb,
                Aabb::new_invalid(),
                Aabb::new_invalid(),
                Aabb::new_invalid(),
            ]),
            children: [id, u32::MAX, u32::MAX, u32::MAX],
            parent: NodeIndex::invalid(),
            flags: QbvhNodeFlags::default(),
        };
    }

    fn do_recurse_rebalance(
        &mut self,
        indices: &mut [usize],
        workspace: &QbvhUpdateWorkspace,
        parent: NodeIndex,
        margin: Real,
    ) -> (u32, Aabb) {
        if indices.len() <= 4 {
            // Leaf case.
            let mut has_leaf = false;
            let mut has_internal = false;

            for id in &*indices {
                has_leaf |= workspace.is_leaf[*id];
                has_internal |= !workspace.is_leaf[*id];
            }

            let my_internal_id = if has_internal {
                self.free_list.pop().unwrap_or_else(|| {
                    self.nodes.push(QbvhNode::empty());
                    self.nodes.len() as u32 - 1
                })
            } else {
                u32::MAX
            };

            let my_leaf_id = if has_leaf {
                self.free_list.pop().unwrap_or_else(|| {
                    self.nodes.push(QbvhNode::empty());
                    self.nodes.len() as u32 - 1
                })
            } else {
                u32::MAX
            };

            let mut my_leaf_aabb = Aabb::new_invalid();
            let mut my_internal_aabb = Aabb::new_invalid();

            let mut leaf_aabbs = [Aabb::new_invalid(); 4];
            let mut internal_aabbs = [Aabb::new_invalid(); 4];

            let mut proxy_ids = [u32::MAX; 4];
            let mut internal_ids = [u32::MAX; 4];
            let mut new_internal_lane_containing_leaf = usize::MAX;

            for (k, id) in indices.iter().enumerate() {
                let orig_id = workspace.orig_ids[*id];
                if workspace.is_leaf[*id] {
                    new_internal_lane_containing_leaf = k;
                    my_leaf_aabb.merge(&workspace.aabbs[*id]);
                    leaf_aabbs[k] = workspace.aabbs[*id];
                    proxy_ids[k] = orig_id;
                    self.proxies[orig_id as usize].node = NodeIndex::new(my_leaf_id, k as u8);
                } else {
                    my_internal_aabb.merge(&workspace.aabbs[*id]);
                    internal_aabbs[k] = workspace.aabbs[*id];
                    internal_ids[k] = orig_id;
                    self.nodes[orig_id as usize].parent = NodeIndex::new(my_internal_id, k as u8);
                }
            }

            if has_internal {
                if has_leaf {
                    internal_ids[new_internal_lane_containing_leaf] = my_leaf_id;
                    internal_aabbs[new_internal_lane_containing_leaf] = my_leaf_aabb;
                    // need to merge the leaf aabb into the internal aabb which we did not do previously
                    my_internal_aabb.merge(&my_leaf_aabb);
                }

                let internal_node = QbvhNode {
                    simd_aabb: SimdAabb::from(internal_aabbs),
                    children: internal_ids,
                    parent,
                    flags: QbvhNodeFlags::default(),
                };
                // NOTE: we don’t dilate because the AABBs for the rebalancing were
                //       already dilated during the refit.
                self.nodes[my_internal_id as usize] = internal_node;
            }

            if has_leaf {
                let leaf_node = QbvhNode {
                    simd_aabb: SimdAabb::from(leaf_aabbs),
                    children: proxy_ids,
                    parent: if has_internal {
                        NodeIndex::new(my_internal_id, new_internal_lane_containing_leaf as u8)
                    } else {
                        parent
                    },
                    flags: QbvhNodeFlags::LEAF,
                };
                // NOTE: we don’t dilate because the AABBs for the rebalancing were
                //       already dilated during the refit.
                self.nodes[my_leaf_id as usize] = leaf_node;
            }

            if has_internal {
                return (my_internal_id, my_internal_aabb);
            } else {
                return (my_leaf_id, my_leaf_aabb);
            }
        }

        // Compute the center and variance along each dimension.
        // In 3D we compute the variance to not-subdivide the dimension with lowest variance.
        // Therefore variance computation is not needed in 2D because we only have 2 dimension
        // to split in the first place.
        let mut center = Point::origin();
        #[cfg(feature = "dim3")]
        let mut variance = Vector::zeros();

        let center_denom = 1.0 / (indices.len() as Real);

        for i in &*indices {
            let coords = workspace.aabbs[*i].center().coords;
            center += coords * center_denom;
        }

        #[cfg(feature = "dim3")]
        {
            let variance_denom = 1.0 / ((indices.len() - 1) as Real);
            for i in &*indices {
                let dir_to_center = workspace.aabbs[*i].center() - center;
                variance += dir_to_center.component_mul(&dir_to_center) * variance_denom;
            }
        }

        // Find the axis with minimum variance. This is the axis along
        // which we are **not** subdividing our set.
        #[allow(unused_mut)] // Does not need to be mutable in 2D.
        let mut subdiv_dims = [0, 1];
        #[cfg(feature = "dim3")]
        {
            let min = variance.imin();
            subdiv_dims[0] = (min + 1) % 3;
            subdiv_dims[1] = (min + 2) % 3;
        }

        let node = QbvhNode {
            simd_aabb: SimdAabb::new_invalid(),
            children: [0; 4], // Will be set after the recursive call
            parent,
            flags: QbvhNodeFlags::default(),
        };

        let nid = if let Some(nid) = self.free_list.pop() {
            self.nodes[nid as usize] = node;
            nid
        } else {
            let nid = self.nodes.len();
            self.nodes.push(node);
            nid as u32
        };

        // Split the set along the two subdiv_dims dimensions.
        let splitter = CenterDataSplitter::default();
        let splits =
            splitter.split_dataset_wo_workspace(subdiv_dims, center, indices, &workspace.aabbs);
        let n = [
            NodeIndex::new(nid, 0),
            NodeIndex::new(nid, 1),
            NodeIndex::new(nid, 2),
            NodeIndex::new(nid, 3),
        ];

        let children = [
            self.do_recurse_rebalance(splits[0], workspace, n[0], margin),
            self.do_recurse_rebalance(splits[1], workspace, n[1], margin),
            self.do_recurse_rebalance(splits[2], workspace, n[2], margin),
            self.do_recurse_rebalance(splits[3], workspace, n[3], margin),
        ];

        // Now we know the indices of the child nodes.
        self.nodes[nid as usize].children =
            [children[0].0, children[1].0, children[2].0, children[3].0];
        self.nodes[nid as usize].simd_aabb =
            SimdAabb::from([children[0].1, children[1].1, children[2].1, children[3].1]);
        self.nodes[nid as usize]
            .simd_aabb
            .loosen(SimdReal::splat(margin));

        let my_aabb = self.nodes[nid as usize].simd_aabb.to_merged_aabb();
        (nid, my_aabb)
    }
}
