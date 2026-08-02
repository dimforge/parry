use super::bvh_tree::{BvhNodeIndex, BvhNodeVec, BvhNodeWide};
use super::{Bvh, BvhNode, BvhWorkspace};
use crate::utils::VecMap;
use alloc::vec::Vec;

/// Raw pointers to the refit buffers, shared across parallel refit tasks.
///
/// Safety: tasks write disjoint index ranges (see `refit_recurse_parallel`).
#[cfg(feature = "parallel")]
#[derive(Copy, Clone)]
struct RefitPtrs {
    target: *mut BvhNodeVec,
    leaf_data: *mut VecMap<BvhNodeIndex>,
    parents: *mut Vec<BvhNodeIndex>,
}

#[cfg(feature = "parallel")]
unsafe impl Send for RefitPtrs {}
#[cfg(feature = "parallel")]
unsafe impl Sync for RefitPtrs {}

impl Bvh {
    /// Updates the BVH's internal node AABBs after leaf changes.
    ///
    /// Refitting ensures that every internal node's AABB tightly encloses the AABBs of its
    /// children. This operation is essential after updating leaf positions with
    /// [`insert_or_update_partially`] and is much faster than rebuilding the entire tree.
    ///
    /// In addition to updating AABBs, this method:
    /// - Reorders nodes in depth-first order for better cache locality during queries
    /// - Ensures leaf counts on each node are correct
    /// - Propagates change flags from leaves to ancestors (for change detection)
    ///
    /// # When to Use
    ///
    /// Call `refit` after:
    /// - Bulk updates with [`insert_or_update_partially`]
    /// - Any operation that modifies leaf AABBs without updating ancestor nodes
    /// - When you want to optimize tree layout for better query performance
    ///
    /// **Don't call `refit` after**:
    /// - Regular [`insert`] calls (they already update ancestors)
    /// - [`remove`] calls (they already maintain tree validity)
    ///
    /// # Arguments
    ///
    /// * `workspace` - A reusable workspace to avoid allocations. Can be shared across
    ///   multiple BVH operations for better performance.
    ///
    /// # Performance
    ///
    /// - **Time**: O(n) where n is the number of nodes
    /// - **Space**: O(n) temporary storage in workspace
    /// - Much faster than rebuilding the tree from scratch
    /// - Essential for maintaining good query performance in dynamic scenes
    ///
    /// # Examples
    ///
    /// ## After bulk updates
    ///
    /// ```
    /// # #[cfg(all(feature = "dim3", feature = "f32"))] {
    /// use parry3d::partitioning::{Bvh, BvhWorkspace};
    /// use parry3d::bounding_volume::Aabb;
    /// use parry3d::math::Vector;
    ///
    /// let mut bvh = Bvh::new();
    /// let mut workspace = BvhWorkspace::default();
    ///
    /// // Insert initial objects
    /// for i in 0..100 {
    ///     let aabb = Aabb::new(
    ///         Vector::new(i as f32, 0.0, 0.0),
    ///         Vector::new(i as f32 + 1.0, 1.0, 1.0)
    ///     );
    ///     bvh.insert(aabb, i);
    /// }
    ///
    /// // Update all objects without tree propagation (faster)
    /// for i in 0..100 {
    ///     let offset = 0.1;
    ///     let aabb = Aabb::new(
    ///         Vector::new(i as f32 + offset, 0.0, 0.0),
    ///         Vector::new(i as f32 + 1.0 + offset, 1.0, 1.0)
    ///     );
    ///     bvh.insert_or_update_partially(aabb, i, 0.0);
    /// }
    ///
    /// // Now update the tree in one efficient pass
    /// bvh.refit(&mut workspace);
    /// # }
    /// ```
    ///
    /// ## In a game loop
    ///
    /// ```
    /// # #[cfg(all(feature = "dim3", feature = "f32"))] {
    /// use parry3d::partitioning::{Bvh, BvhWorkspace};
    /// use parry3d::bounding_volume::Aabb;
    /// use parry3d::math::Vector;
    ///
    /// let mut bvh = Bvh::new();
    /// let mut workspace = BvhWorkspace::default();
    ///
    /// // Game initialization - add objects
    /// for i in 0..1000 {
    ///     let aabb = Aabb::new(
    ///         Vector::new(i as f32, 0.0, 0.0),
    ///         Vector::new(i as f32 + 1.0, 1.0, 1.0)
    ///     );
    ///     bvh.insert(aabb, i);
    /// }
    ///
    /// // Game loop - update objects each frame
    /// for frame in 0..100 {
    ///     // Update physics, AI, etc.
    ///     for i in 0..1000 {
    ///         let time = frame as f32 * 0.016; // ~60 FPS
    ///         let pos = time.sin() * 10.0;
    ///         let aabb = Aabb::new(
    ///             Vector::new(i as f32 + pos, 0.0, 0.0),
    ///             Vector::new(i as f32 + pos + 1.0, 1.0, 1.0)
    ///         );
    ///         bvh.insert_or_update_partially(aabb, i, 0.0);
    ///     }
    ///
    ///     // Refit once per frame for all updates
    ///     bvh.refit(&mut workspace);
    ///
    ///     // Now perform collision detection queries...
    /// }
    /// # }
    /// ```
    ///
    /// ## With change detection margin
    ///
    /// ```
    /// # #[cfg(all(feature = "dim3", feature = "f32"))] {
    /// use parry3d::partitioning::{Bvh, BvhWorkspace};
    /// use parry3d::bounding_volume::Aabb;
    /// use parry3d::math::Vector;
    ///
    /// let mut bvh = Bvh::new();
    /// let mut workspace = BvhWorkspace::default();
    ///
    /// // Add an object
    /// let aabb = Aabb::new(Vector::ZERO, Vector::new(1.0, 1.0, 1.0));
    /// bvh.insert(aabb, 0);
    ///
    /// // Update with a margin - tree won't update if movement is small
    /// let margin = 0.5;
    /// let new_aabb = Aabb::new(Vector::new(0.1, 0.0, 0.0), Vector::new(1.1, 1.0, 1.0));
    /// bvh.insert_or_update_partially(new_aabb, 0, margin);
    ///
    /// // Refit propagates the change detection flags
    /// bvh.refit(&mut workspace);
    /// # }
    /// ```
    ///
    /// # Comparison with `refit_without_opt`
    ///
    /// This method reorganizes the tree in memory for better cache performance.
    /// If you only need to update AABBs without reordering, use [`refit_without_opt`](Self::refit_without_opt)
    /// which is faster but doesn't improve memory layout.
    ///
    /// # Notes
    ///
    /// - Reuses the provided `workspace` to avoid allocations
    /// - Safe to call even if no leaves were modified (just reorganizes tree)
    /// - Does not change the tree's topology, only AABBs and layout
    /// - Call this before [`optimize_incremental`] for best results
    ///
    /// # See Also
    ///
    /// - [`insert_or_update_partially`](Bvh::insert_or_update_partially) - Update leaves
    ///   without propagation
    /// - [`refit_without_opt`](Self::refit_without_opt) - Faster refit without memory
    ///   reorganization
    /// - [`optimize_incremental`](Bvh::optimize_incremental) - Improve tree quality
    /// - [`BvhWorkspace`] - Reusable workspace for operations
    ///
    /// [`insert_or_update_partially`]: Bvh::insert_or_update_partially
    /// [`insert`]: Bvh::insert
    /// [`remove`]: Bvh::remove
    /// [`optimize_incremental`]: Bvh::optimize_incremental
    pub fn refit(&mut self, workspace: &mut BvhWorkspace) {
        Self::refit_buffers::<true>(
            &mut self.nodes,
            &mut workspace.refit_tmp,
            &mut self.leaf_node_indices,
            &mut self.parents,
        );

        // Swap the old nodes with the refitted ones.
        core::mem::swap(&mut self.nodes, &mut workspace.refit_tmp);
        // The refit rebuilt the node array in depth-first order, dropping the
        // orphaned slots the free list pointed to.
        self.free_wide_nodes.clear();
    }

    /// Same as [`Self::refit`], but processes independent subtrees in parallel.
    ///
    /// The result is identical to [`Self::refit`] (same node layout, same flags);
    /// only the work distribution differs.
    #[cfg(feature = "parallel")]
    pub fn refit_parallel(&mut self, workspace: &mut BvhWorkspace) {
        Self::refit_buffers_parallel::<true>(
            &mut self.nodes,
            &mut workspace.refit_tmp,
            &mut self.leaf_node_indices,
            &mut self.parents,
        );

        // Swap the old nodes with the refitted ones.
        core::mem::swap(&mut self.nodes, &mut workspace.refit_tmp);
        // The refit rebuilt the node array in depth-first order, dropping the
        // orphaned slots the free list pointed to.
        self.free_wide_nodes.clear();
    }

    #[cfg(feature = "parallel")]
    fn refit_buffers_parallel<const RESOLVE: bool>(
        source: &mut BvhNodeVec,
        target: &mut BvhNodeVec,
        leaf_data: &mut VecMap<BvhNodeIndex>,
        parents: &mut Vec<BvhNodeIndex>,
    ) {
        // Subtrees below this leaf count are refitted sequentially.
        const SEQ_LEAF_THRESHOLD: u32 = 2048;
        // Bounds the number of spawned tasks to 2^MAX_SPLIT_DEPTH.
        const MAX_SPLIT_DEPTH: u32 = 6;

        if source.is_empty() || source[0].leaf_count() <= SEQ_LEAF_THRESHOLD.max(2) {
            return Self::refit_buffers::<RESOLVE>(source, target, leaf_data, parents);
        }

        target.resize(
            source.len(),
            BvhNodeWide {
                left: BvhNode::zeros(),
                right: BvhNode::zeros(),
            },
        );
        parents.resize(source.len(), BvhNodeIndex::default());

        let ptrs = RefitPtrs {
            target: target as *mut BvhNodeVec,
            leaf_data: leaf_data as *mut VecMap<BvhNodeIndex>,
            parents: parents as *mut Vec<BvhNodeIndex>,
        };

        // Mirror of `refit_buffers`' root special case, with the sequential target-id
        // counter replaced by offsets computed from the subtree leaf counts: the
        // depth-first layout of a subtree with `n` leaves spans exactly `n - 1` wide
        // nodes.
        let root = source[0];
        let left_size = if root.left.is_leaf() {
            0
        } else {
            root.left.leaf_count() - 1
        };
        let right_size = if root.right.is_leaf() {
            0
        } else {
            root.right.leaf_count() - 1
        };
        let final_len = (1 + left_size + right_size) as usize;

        let _ = rayon::join(
            || {
                // SAFETY: this task writes only to target slots [1, 1 + left_size),
                //         entry target[0].left, and the leaf data of leaves of the
                //         left subtree — all disjoint from the other task.
                let target = unsafe { &mut *ptrs.target };
                if !root.left.is_leaf() {
                    Self::refit_recurse_parallel::<RESOLVE>(
                        source,
                        ptrs,
                        root.left.children,
                        1,
                        BvhNodeIndex::left(0),
                        MAX_SPLIT_DEPTH,
                        SEQ_LEAF_THRESHOLD,
                    );
                } else {
                    target[0].left = root.left;
                    if RESOLVE {
                        target[0].left.data.resolve_pending_change();
                    }
                }
            },
            || {
                // SAFETY: see the other task; slots [1 + left_size, final_len) and
                //         entry target[0].right.
                let target = unsafe { &mut *ptrs.target };
                if !root.right.is_leaf() {
                    Self::refit_recurse_parallel::<RESOLVE>(
                        source,
                        ptrs,
                        root.right.children,
                        1 + left_size,
                        BvhNodeIndex::right(0),
                        MAX_SPLIT_DEPTH,
                        SEQ_LEAF_THRESHOLD,
                    );
                } else {
                    target[0].right = root.right;
                    if RESOLVE {
                        target[0].right.data.resolve_pending_change();
                    }
                }
            },
        );

        source.truncate(final_len);
        target.truncate(final_len);
        parents.truncate(final_len);
    }

    /// Recursive parallel counterpart of `refit_recurse`.
    ///
    /// `target_id` is the depth-first slot this subtree's root occupies (its left
    /// child subtree starts at `target_id + 1`, its right child subtree right after
    /// the left one, whose extent is known from its leaf count).
    #[cfg(feature = "parallel")]
    fn refit_recurse_parallel<const RESOLVE: bool>(
        source: &BvhNodeVec,
        ptrs: RefitPtrs,
        source_id: u32,
        target_id: u32,
        parent: BvhNodeIndex,
        depth: u32,
        seq_leaf_threshold: u32,
    ) {
        let node = source[source_id as usize];
        let leaf_count = node.left.leaf_count() + node.right.leaf_count();

        if depth == 0 || leaf_count <= seq_leaf_threshold {
            //

            // SAFETY: the sequential refit of this subtree only touches the target
            //         slots [target_id, target_id + leaf_count - 1), its parent entry,
            //         and its own leaves' data: all disjoint from the other tasks.
            let target = unsafe { &mut *ptrs.target };
            let leaf_data = unsafe { &mut *ptrs.leaf_data };
            let parents = unsafe { &mut *ptrs.parents };
            let mut counter = target_id;
            Self::refit_recurse::<RESOLVE>(
                source,
                target,
                leaf_data,
                parents,
                source_id,
                &mut counter,
                parent,
            );
            debug_assert_eq!(counter, target_id + leaf_count - 1);
            return;
        }

        let left_size = if node.left.is_leaf() {
            0
        } else {
            node.left.leaf_count() - 1
        };

        let _ = rayon::join(
            || {
                let target = unsafe { &mut *ptrs.target };
                let leaf_data = unsafe { &mut *ptrs.leaf_data };
                if !node.left.is_leaf() {
                    Self::refit_recurse_parallel::<RESOLVE>(
                        source,
                        ptrs,
                        node.left.children,
                        target_id + 1,
                        BvhNodeIndex::left(target_id),
                        depth - 1,
                        seq_leaf_threshold,
                    );
                } else {
                    target[target_id as usize].left = node.left;
                    if RESOLVE {
                        target[target_id as usize]
                            .left
                            .data
                            .resolve_pending_change();
                    }
                    leaf_data[node.left.children as usize] = BvhNodeIndex::left(target_id);
                }
            },
            || {
                let target = unsafe { &mut *ptrs.target };
                let leaf_data = unsafe { &mut *ptrs.leaf_data };
                if !node.right.is_leaf() {
                    Self::refit_recurse_parallel::<RESOLVE>(
                        source,
                        ptrs,
                        node.right.children,
                        target_id + 1 + left_size,
                        BvhNodeIndex::right(target_id),
                        depth - 1,
                        seq_leaf_threshold,
                    );
                } else {
                    target[target_id as usize].right = node.right;
                    if RESOLVE {
                        target[target_id as usize]
                            .right
                            .data
                            .resolve_pending_change();
                    }
                    leaf_data[node.right.children as usize] = BvhNodeIndex::right(target_id);
                }
            },
        );

        // Both children of this wide node are now written: compute the summary entry
        // in the parent, like the tail of `refit_recurse`.
        let target = unsafe { &mut *ptrs.target };
        let parents = unsafe { &mut *ptrs.parents };
        let merged = target[target_id as usize]
            .left
            .merged(&target[target_id as usize].right, target_id);
        target[parent] = merged;
        parents[target_id as usize] = parent;
    }

    /// Same as [`Self::refit`], but leaves every change-detection flag untouched.
    ///
    /// Use this to make the tree valid for queries after a batch of
    /// [`Self::insert_or_update_partially`] without consuming the pending change
    /// flags: a later flag-resolving [`Self::refit`] (or
    /// [`Self::refit_partial`]) will promote them as if this call never happened.
    pub fn refit_without_resolve(&mut self, workspace: &mut BvhWorkspace) {
        Self::refit_buffers::<false>(
            &mut self.nodes,
            &mut workspace.refit_tmp,
            &mut self.leaf_node_indices,
            &mut self.parents,
        );
        core::mem::swap(&mut self.nodes, &mut workspace.refit_tmp);
        // The refit rebuilt the node array in depth-first order, dropping the
        // orphaned slots the free list pointed to.
        self.free_wide_nodes.clear();
    }

    /// Parallel version of [`Self::refit_without_resolve`].
    #[cfg(feature = "parallel")]
    pub fn refit_without_resolve_parallel(&mut self, workspace: &mut BvhWorkspace) {
        Self::refit_buffers_parallel::<false>(
            &mut self.nodes,
            &mut workspace.refit_tmp,
            &mut self.leaf_node_indices,
            &mut self.parents,
        );
        core::mem::swap(&mut self.nodes, &mut workspace.refit_tmp);
        // The refit rebuilt the node array in depth-first order, dropping the
        // orphaned slots the free list pointed to.
        self.free_wide_nodes.clear();
    }

    pub(super) fn refit_buffers<const RESOLVE: bool>(
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
            if RESOLVE {
                target[0].left.data.resolve_pending_change();
                if target[0].right.leaf_count() > 0 {
                    target[0].right.data.resolve_pending_change();
                }
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
                Self::refit_recurse::<RESOLVE>(
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
                if RESOLVE {
                    target[0].left.data.resolve_pending_change();
                }

                // NOTE: updating the leaf_data shouldn’t be needed here since the root
                //       is always at 0.
                // *self.leaf_data.get_mut_unknown_gen(left_child_id).unwrap() = BvhNodeIndex::left(0);
            }

            if !source[0].right.is_leaf() {
                Self::refit_recurse::<RESOLVE>(
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
                if RESOLVE {
                    target[0].right.data.resolve_pending_change();
                }
                // NOTE: updating the leaf_data shouldn’t be needed here since the root
                //       is always at 0.
                // *self.leaf_data.get_mut_unknown_gen(right_child_id).unwrap() = BvhNodeIndex::right(0);
            }

            source.truncate(len as usize);
            target.truncate(len as usize);
            parents.truncate(len as usize);
        }
    }

    fn refit_recurse<const RESOLVE: bool>(
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
            Self::refit_recurse::<RESOLVE>(
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
            if RESOLVE {
                target[target_id as usize]
                    .left
                    .data
                    .resolve_pending_change();
            }
            leaf_data[node.left.children as usize] = BvhNodeIndex::left(target_id);
        }

        if !right_is_leaf {
            Self::refit_recurse::<RESOLVE>(
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
            if RESOLVE {
                target[target_id as usize]
                    .right
                    .data
                    .resolve_pending_change();
            }
            leaf_data[node.right.children as usize] = BvhNodeIndex::right(target_id);
        }

        let node = &target[target_id as usize];
        target[parent] = node.left.merged(&node.right, target_id);
        parents[target_id as usize] = parent;
    }

    /// Incrementally refits the tree after a small number of leaf updates or
    /// insertions.
    ///
    /// This is a faster alternative to [`Self::refit`] valid only if, since the last
    /// refit, the only tree modifications were calls to
    /// [`Self::insert_or_update_partially`], [`Self::insert`],
    /// [`Self::insert_with_change_detection`], or
    /// [`Self::reinsert_or_update_with_change_detection`] (in-place updates,
    /// insertions, or removal-based relocations of leaves). In particular, no leaf
    /// was removed without being re-inserted in the same batch, and no optimization
    /// ran. Otherwise, call [`Self::refit`] instead.
    ///
    /// Insertions are safe here because `insert_new_unchecked` keeps the tree
    /// geometrically valid on its own (it enlarges the ancestor AABBs and increments
    /// their leaf counts during its descent) and creates the new leaf with a pending
    /// change flag. Any transient change-flag state it leaves behind (the pending
    /// flag inherited by the wide node that used to hold the insertion sibling, the
    /// raw-merged flags written by its SAH rotations) lies on the inserted leaf's
    /// ancestor path, which the walk below rewrites all the way to the root — see
    /// `refit_path`.
    ///
    /// [`BvhLeafUpdateStatus::Inserted`]: super::BvhLeafUpdateStatus::Inserted
    ///
    /// `previously_changed` must contain (a superset of) the leaves whose change flag
    /// was set by the previous refit: their change flag gets cleared. `newly_changed`
    /// must contain (a superset of) the leaves updated in-place or inserted since the
    /// last refit: their change flag gets set if their fat AABB actually changed
    /// (inserted leaves always count as changed).
    ///
    /// Unlike [`Self::refit`], this runs in `O(changed * tree_height)` instead of
    /// `O(node_count)`, but doesn't reorder nodes in memory.
    pub fn refit_partial(&mut self, previously_changed: &[u32], newly_changed: &[u32]) {
        // First resolve the change flags of every impacted leaf, and only then walk
        // their ancestor paths. Walking while some leaves still hold an unresolved
        // pending flag would propagate that transient state into internal nodes
        // (a pending internal node reads as "unchanged" during traversals).
        for leaf in previously_changed {
            let Some(leaf_node_id) = self.leaf_node_indices.get(*leaf as usize).copied() else {
                continue;
            };
            let data = &mut self.nodes[leaf_node_id].data;

            // Clear the flag of leaves that were changed at the previous refit and
            // didn't move since. Leaves that moved again (pending) are promoted by
            // the next loop instead.
            if !data.is_change_pending() && data.is_changed() {
                data.resolve_pending_change();
            }
        }

        for leaf in newly_changed {
            let Some(leaf_node_id) = self.leaf_node_indices.get(*leaf as usize).copied() else {
                continue;
            };
            let data = &mut self.nodes[leaf_node_id].data;

            // Promote pending leaves to CHANGED. Leaves whose update stayed within
            // the change-detection margin have no flag to update.
            if data.is_change_pending() {
                data.resolve_pending_change();
            }
        }

        for leaf in previously_changed.iter().chain(newly_changed) {
            let Some(leaf_node_id) = self.leaf_node_indices.get(*leaf as usize).copied() else {
                continue;
            };
            self.refit_path(leaf_node_id);
        }
    }

    /// Propagates AABB and change-flag updates from the given node up to the root.
    ///
    /// The walk deliberately does NOT stop early when an ancestor's recomputed value
    /// matches its stored one: leaf insertions can leave transient change-flag state
    /// (pending flags inherited by the wide node that used to hold the insertion
    /// sibling, raw-merged flags written by the insertion's SAH rotations) anywhere
    /// on the inserted leaf's ancestor path. Walking the whole path rewrites every
    /// ancestor with an exact, normalized recomputation, which is what keeps the
    /// internal change flags exactly equal to the OR of their descendant leaves'
    /// flags. An early-out could strand a stale pending flag above the break point,
    /// and a pending internal node reads as "unchanged" during change-detection
    /// traversals — hiding every changed leaf underneath. The full walk costs
    /// O(tree height) per changed leaf, which is fine in the small-change regime
    /// partial refits are meant for.
    fn refit_path(&mut self, node: BvhNodeIndex) {
        let (mut wide_id, _) = node.decompose();

        while wide_id != 0 {
            let parent = self.parents[wide_id];
            let wide = &self.nodes[wide_id];
            let mut recomputed = wide.left.merged(&wide.right, wide_id as u32);
            recomputed.data.normalize_change_flag();

            self.nodes[parent] = recomputed;
            (wide_id, _) = parent.decompose();
        }
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
