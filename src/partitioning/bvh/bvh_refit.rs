use super::bvh_tree::{BvhNodeIndex, BvhNodeVec, BvhNodeWide};
use super::{Bvh, BvhNode, BvhWorkspace};
use alloc::vec::Vec;
use vec_map::VecMap;

impl Bvh {
    /// Performs a tree refitting with internal storage cache optimizations.
    ///
    /// Refitting is the action of ensuring that every node’s AABB encloses the AABBs of its
    /// children. In addition, this method will:
    /// - Ensure that nodes are stored internally in depth-first order for better cache locality
    ///   on depth-first searches.
    /// - Ensure that the leaf count on each node is correct.
    /// - Propagate the [`BvhNode::is_changed`] flag changes detected by [`Self::insert_or_update_partially`]
    ///   (if `change_detection_margin` was nonzero) from leaves to its ascendants.
    pub fn refit(&mut self, workspace: &mut BvhWorkspace) {
        Self::refit_buffers(
            &mut self.nodes,
            &mut workspace.refit_tmp,
            &mut self.leaf_node_indices,
            &mut self.parents,
        );

        // Swap the old nodes with the refitted ones.
        core::mem::swap(&mut self.nodes, &mut workspace.refit_tmp);
    }

    pub(super) fn refit_buffers(
        source: &mut BvhNodeVec,
        target: &mut BvhNodeVec,
        leaf_data: &mut VecMap<BvhNodeIndex>,
        parents: &mut Vec<BvhNodeIndex>,
    ) {
        if source.is_empty() {
            target.clear();
            parents.clear();
        } else if source[0].leaf_count() <= 2 {
            // No actual refit to apply, just copy the root wide node.
            target.clear();
            parents.clear();
            target.push(source[0]);
            target[0].left.data.resolve_pending_change();
            if target[0].right.leaf_count() > 0 {
                target[0].right.data.resolve_pending_change();
            }
            parents.push(BvhNodeIndex::default());
        } else if !source.is_empty() && source[0].leaf_count() > 2 {
            target.resize(
                source.len(),
                BvhNodeWide {
                    left: BvhNode::zeros(),
                    right: BvhNode::zeros(),
                },
            );

            let mut len = 1;

            // Start with a special case for the root then recurse.
            let left_child_id = source[0].left.children;
            let right_child_id = source[0].right.children;

            if !source[0].left.is_leaf() {
                Self::refit_recurse(
                    source,
                    target,
                    leaf_data,
                    parents,
                    left_child_id,
                    &mut len,
                    BvhNodeIndex::left(0),
                );
            } else {
                target[0].left = source[0].left;
                target[0].left.data.resolve_pending_change();

                // NOTE: updating the leaf_data shouldn’t be needed here since the root
                //       is always at 0.
                // *self.leaf_data.get_mut_unknown_gen(left_child_id).unwrap() = BvhNodeIndex::left(0);
            }

            if !source[0].right.is_leaf() {
                Self::refit_recurse(
                    source,
                    target,
                    leaf_data,
                    parents,
                    right_child_id,
                    &mut len,
                    BvhNodeIndex::right(0),
                );
            } else {
                target[0].right = source[0].right;
                target[0].right.data.resolve_pending_change();
                // NOTE: updating the leaf_data shouldn’t be needed here since the root
                //       is always at 0.
                // *self.leaf_data.get_mut_unknown_gen(right_child_id).unwrap() = BvhNodeIndex::right(0);
            }

            source.truncate(len as usize);
            target.truncate(len as usize);
            parents.truncate(len as usize);
        }
    }

    fn refit_recurse(
        source: &BvhNodeVec,
        target: &mut BvhNodeVec,
        leaf_data: &mut VecMap<BvhNodeIndex>,
        parents: &mut [BvhNodeIndex],
        source_id: u32,
        target_id_mut: &mut u32,
        parent: BvhNodeIndex,
    ) {
        let target_id = *target_id_mut;
        *target_id_mut += 1;

        let node = &source[source_id as usize];
        let left_is_leaf = node.left.is_leaf();
        let right_is_leaf = node.right.is_leaf();
        let left_source_id = node.left.children;
        let right_source_id = node.right.children;

        if !left_is_leaf {
            Self::refit_recurse(
                source,
                target,
                leaf_data,
                parents,
                left_source_id,
                target_id_mut,
                BvhNodeIndex::left(target_id),
            );
        } else {
            let node = &source[source_id as usize];
            target[target_id as usize].left = node.left;
            target[target_id as usize]
                .left
                .data
                .resolve_pending_change();
            leaf_data[node.left.children as usize] = BvhNodeIndex::left(target_id);
        }

        if !right_is_leaf {
            Self::refit_recurse(
                source,
                target,
                leaf_data,
                parents,
                right_source_id,
                target_id_mut,
                BvhNodeIndex::right(target_id),
            );
        } else {
            let node = &source[source_id as usize];
            target[target_id as usize].right = node.right;
            target[target_id as usize]
                .right
                .data
                .resolve_pending_change();
            leaf_data[node.right.children as usize] = BvhNodeIndex::right(target_id);
        }

        let node = &target[target_id as usize];
        target[parent] = node.left.merged(&node.right, target_id);
        parents[target_id as usize] = parent;
    }

    /// Similar to [`Self::refit`] but without any optimization of the internal node storage layout.
    ///
    /// This can be faster than [`Self::refit`] but doesn’t reorder node to be more cache-efficient
    /// on tree traversals.
    pub fn refit_without_opt(&mut self) {
        if self.leaf_count() > 2 {
            let root = &self.nodes[0];
            let left = root.left.children;
            let right = root.right.children;
            let left_is_leaf = root.left.is_leaf();
            let right_is_leaf = root.right.is_leaf();

            if !left_is_leaf {
                self.recurse_refit_without_opt(left, BvhNodeIndex::left(0));
            }
            if !right_is_leaf {
                self.recurse_refit_without_opt(right, BvhNodeIndex::right(0));
            }
        }
    }

    fn recurse_refit_without_opt(&mut self, node_id: u32, parent: BvhNodeIndex) {
        let node = &self.nodes[node_id as usize];
        let left = &node.left;
        let right = &node.right;
        let left_is_leaf = left.is_leaf();
        let right_is_leaf = right.is_leaf();
        let left_children = left.children;
        let right_children = right.children;

        if !left_is_leaf {
            self.recurse_refit_without_opt(left_children, BvhNodeIndex::left(node_id));
        } else {
            self.nodes[node_id as usize]
                .left
                .data
                .resolve_pending_change();
        }
        if !right_is_leaf {
            self.recurse_refit_without_opt(right_children, BvhNodeIndex::right(node_id));
        } else {
            self.nodes[node_id as usize]
                .right
                .data
                .resolve_pending_change();
        }

        let node = &self.nodes[node_id as usize];
        let left = &node.left;
        let right = &node.right;
        let merged = left.merged(right, node_id);

        self.nodes[parent] = merged;
    }
}
