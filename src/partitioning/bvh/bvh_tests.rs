use crate::bounding_volume::Aabb;
use crate::math::{Real, Vector};
use crate::partitioning::{
    Bvh, BvhBuildStrategy, BvhNode, BvhNodeIndex, BvhWorkspace, TraversalAction,
};

fn make_test_aabb(i: usize) -> Aabb {
    Aabb::from_half_extents(Vector::splat(i as Real).into(), Vector::splat(1.0))
}

#[test]
fn test_leaves_iteration() {
    let leaves = [
        make_test_aabb(0), // mins at (0,0,0) - should pass
        make_test_aabb(5), // mins at (5,5,5) - should be filtered out
    ];
    let bvh = Bvh::from_leaves(BvhBuildStrategy::Binned, &leaves);

    // Only allow nodes with mins.x <= 3.0 (should only pass leaf 0)
    let check = |node: &BvhNode| -> bool { node.mins.x <= 3.0 };

    let mut found_invalid_leaf = false;
    for leaf_index in bvh.leaves(check) {
        if leaf_index == 1 {
            // This is the leaf that should be filtered out
            found_invalid_leaf = true;
            break;
        }
    }

    if found_invalid_leaf {
        panic!("Leaves iterator returned an invalid leaf");
    }
}

#[test]
fn test_traverse_indexed() {
    // Empty tree: callback must never fire, regardless of `subtree` being `None`.
    let empty = Bvh::new();
    empty.traverse_indexed(None, |_, _| {
        panic!("callback should not be called on an empty BVH");
    });

    // Single-leaf tree exercises the partial-root branch when starting from the root.
    let single = Bvh::from_leaves(BvhBuildStrategy::Binned, &[make_test_aabb(0)]);
    let mut single_visited = std::vec::Vec::new();
    single.traverse_indexed(None, |node, idx| {
        single_visited.push((idx, node.leaf_data()));
        TraversalAction::Continue
    });
    assert_eq!(single_visited.len(), 1);
    assert_eq!(single_visited[0].0, BvhNodeIndex::left(0));
    assert_eq!(single_visited[0].1, Some(0));

    // Multi-leaf tree: traversing from the root must visit every leaf, and the
    // index passed to the callback must round-trip through `bvh.nodes`.
    let leaves: std::vec::Vec<_> = (0..16).map(make_test_aabb).collect();
    let bvh = Bvh::from_leaves(BvhBuildStrategy::Binned, &leaves);

    let mut seen_leaves = std::vec::Vec::new();
    let mut traverse_indexed_calls = std::vec::Vec::new();
    bvh.traverse_indexed(None, |node, idx| {
        // Every reported index must point to the same node we just received.
        let by_idx: &BvhNode = &bvh.nodes[idx];
        assert!(core::ptr::eq(by_idx, node));

        traverse_indexed_calls.push(idx);
        if let Some(data) = node.leaf_data() {
            seen_leaves.push(data);
        }
        TraversalAction::Continue
    });
    seen_leaves.sort();
    assert_eq!(seen_leaves, (0..16).collect::<std::vec::Vec<_>>());

    // `traverse_indexed(None, ...)` must visit exactly the same nodes (in the same
    // order) as `traverse`.
    let mut traverse_nodes: std::vec::Vec<*const BvhNode> = std::vec::Vec::new();
    bvh.traverse(|node| {
        traverse_nodes.push(node as *const _);
        TraversalAction::Continue
    });
    let indexed_nodes: std::vec::Vec<*const BvhNode> = traverse_indexed_calls
        .iter()
        .map(|idx| &bvh.nodes[*idx] as *const _)
        .collect();
    assert_eq!(traverse_nodes, indexed_nodes);

    // Starting from a specific subtree must only visit that subtree (the start
    // node and its descendants), and every reported leaf must belong to it.
    let subtree_root_idx = BvhNodeIndex::left(0);
    let mut subtree_leaves = std::vec::Vec::new();
    let mut subtree_visited = std::vec::Vec::new();
    bvh.traverse_indexed(Some(subtree_root_idx), |node, idx| {
        subtree_visited.push(idx);
        if let Some(data) = node.leaf_data() {
            subtree_leaves.push(data);
        }
        TraversalAction::Continue
    });
    assert_eq!(subtree_visited[0], subtree_root_idx);
    // The subtree's leaves must be a non-empty strict subset of the full set.
    assert!(!subtree_leaves.is_empty());
    assert!(subtree_leaves.len() < 16);
    for leaf in &subtree_leaves {
        assert!(seen_leaves.contains(leaf));
    }
    // Leaf count reported by the subtree's root must match the visited leaves.
    assert_eq!(
        bvh.nodes[subtree_root_idx].leaf_count() as usize,
        subtree_leaves.len()
    );

    // Starting from a leaf node visits exactly that leaf.
    let leaf_idx = *traverse_indexed_calls
        .iter()
        .find(|idx| bvh.nodes[**idx].is_leaf())
        .expect("the tree must contain at least one leaf");
    let mut leaf_only = std::vec::Vec::new();
    bvh.traverse_indexed(Some(leaf_idx), |node, idx| {
        leaf_only.push((idx, node.leaf_data()));
        TraversalAction::Continue
    });
    assert_eq!(leaf_only.len(), 1);
    assert_eq!(leaf_only[0].0, leaf_idx);
    assert!(leaf_only[0].1.is_some());

    // `Prune` at the start node must visit it once and stop.
    let mut prune_visits = 0;
    bvh.traverse_indexed(Some(BvhNodeIndex::left(0)), |_, _| {
        prune_visits += 1;
        TraversalAction::Prune
    });
    assert_eq!(prune_visits, 1);

    // `EarlyExit` at the start node must visit it once and stop.
    let mut exit_visits = 0;
    bvh.traverse_indexed(Some(BvhNodeIndex::left(0)), |_, _| {
        exit_visits += 1;
        TraversalAction::EarlyExit
    });
    assert_eq!(exit_visits, 1);

    // `EarlyExit` partway through must short-circuit the full traversal.
    let mut early = 0;
    bvh.traverse_indexed(None, |_, _| {
        early += 1;
        if early >= 3 {
            TraversalAction::EarlyExit
        } else {
            TraversalAction::Continue
        }
    });
    assert_eq!(early, 3);
}

#[test]
fn bvh_build_and_removal() {
    // Check various combination of building pattern and removal pattern.
    // The tree validity is asserted at every step.
    #[derive(Copy, Clone, Debug)]
    enum BuildPattern {
        Ploc,
        Binned,
        Insert,
    }

    #[derive(Copy, Clone, Debug)]
    enum RemovalPattern {
        InOrder,
        RevOrder,
        EvenOdd,
    }

    for build_pattern in [
        BuildPattern::Ploc,
        BuildPattern::Binned,
        BuildPattern::Insert,
    ] {
        for removal_pattern in [
            RemovalPattern::InOrder,
            RemovalPattern::RevOrder,
            RemovalPattern::EvenOdd,
        ] {
            for len in 1..=100 {
                std::println!(
                    "Testing build: {:?}, removal: {:?}, len: {}",
                    build_pattern,
                    removal_pattern,
                    len
                );
                let leaves: std::vec::Vec<_> = (0..len).map(make_test_aabb).collect();

                let mut bvh = match build_pattern {
                    BuildPattern::Binned => Bvh::from_leaves(BvhBuildStrategy::Binned, &leaves),
                    BuildPattern::Ploc => Bvh::from_leaves(BvhBuildStrategy::Ploc, &leaves),
                    BuildPattern::Insert => {
                        let mut bvh = Bvh::new();
                        for i in 0..len {
                            bvh.insert(make_test_aabb(i), i as u32);
                            bvh.assert_well_formed();
                        }
                        bvh
                    }
                };

                for _ in 0..3 {
                    bvh.assert_well_formed();

                    match removal_pattern {
                        RemovalPattern::InOrder => {
                            // Remove in insertion order.
                            for i in 0..len {
                                bvh.remove(i as u32);
                                bvh.assert_well_formed();
                            }
                        }
                        RemovalPattern::RevOrder => {
                            // Remove in reverse insertion order.
                            for i in (0..len).rev() {
                                bvh.remove(i as u32);
                                bvh.assert_well_formed();
                            }
                        }
                        RemovalPattern::EvenOdd => {
                            // Remove even indices first, then odd.
                            for i in (0..len).filter(|i| i % 2 == 0) {
                                bvh.remove(i as u32);
                                bvh.assert_well_formed();
                            }
                            for i in (0..len).filter(|i| i % 2 != 0) {
                                bvh.remove(i as u32);
                                bvh.assert_well_formed();
                            }
                        }
                    }

                    // Re-insert everything.
                    for (i, leaf) in leaves.iter().enumerate() {
                        bvh.insert(*leaf, i as u32);
                    }
                }
            }
        }
    }
}

#[test]
fn bvh_remove_to_partial_root_then_optimize() {
    // Regression test for the bug reported in #409 where `Bvh::remove` would leave orphaned
    // wide nodes in `self.nodes`/`self.parents` after reducing the tree to a
    // partial root (a single surviving leaf at node 0). Earlier removals on a
    // tree with more than one wide node would intentionally leave orphans for
    // the next `refit` to compact, but if a partial root was created before
    // any refit, those orphans remained reachable from `self.nodes` and
    // `optimize_incremental` would walk them as if they were live, crashing
    // on the corrupt structure.
    //
    // We pick enough leaves to ensure at least one orphan-leaving removal
    // (`wide_node_index != 0`) before the final partial-root removal.
    let leaves: std::vec::Vec<_> = (0..10).map(make_test_aabb).collect();
    let mut bvh = Bvh::from_leaves(BvhBuildStrategy::Binned, &leaves);

    // Remove all but the last leaf, without ever calling refit in between.
    for i in 0..(leaves.len() as u32 - 1) {
        bvh.remove(i);
    }

    // After the final remove, the tree should be a partial root with exactly
    // one surviving leaf, and no orphaned wide nodes left over.
    assert_eq!(bvh.leaf_count(), 1);

    // Without the fix this call walks the orphan-laden tree as if it were
    // live and ends up corrupting/crashing on the partial root.
    let mut workspace = BvhWorkspace::default();
    bvh.optimize_incremental(&mut workspace);

    assert_eq!(bvh.nodes.len(), 1);
    assert_eq!(bvh.parents.len(), 1);
    bvh.assert_well_formed();
}

#[cfg(feature = "parallel")]
#[cfg(all(feature = "dim3", feature = "f32"))]
mod parallel_batch_update {
    use crate::bounding_volume::Aabb;
    use crate::math::Vector;
    use crate::partitioning::{Bvh, BvhLeafUpdateStatus, BvhWorkspace};
    use alloc::vec::Vec;

    /// The parallel batch leaf update must behave exactly like the sequential
    /// per-leaf `insert_or_update_partially` calls (statuses and final tree).
    #[test]
    fn batch_partial_update_matches_sequential() {
        let aabb = |i: usize, offset: f32| {
            Aabb::new(
                Vector::new(i as f32 + offset, 0.0, 0.0).into(),
                Vector::new(i as f32 + offset + 1.0, 1.0, 1.0).into(),
            )
        };

        let mut seq = Bvh::new();
        let mut par = Bvh::new();
        for i in 0..1000 {
            seq.insert(aabb(i, 0.0), i as u32);
            par.insert(aabb(i, 0.0), i as u32);
        }
        // Fatten every leaf so the batch below exercises the within-margin
        // (`Unchanged`) path too.
        for i in 0..1000 {
            let _ = seq.insert_or_update_partially(aabb(i, 2.5), i as u32, 0.5);
            let _ = par.insert_or_update_partially(aabb(i, 2.5), i as u32, 0.5);
        }

        // A mix of: moved-beyond-margin leaves, moved-within-margin leaves, and
        // brand new leaves.
        let updates: Vec<(Aabb, u32, f32)> = (0..1500)
            .map(|i| {
                let offset = if i % 3 == 0 { 2.51 } else { 5.0 };
                (aabb(i, offset), i as u32, 0.1)
            })
            .collect();

        let seq_statuses: Vec<BvhLeafUpdateStatus> = updates
            .iter()
            .map(|(aabb, leaf, margin)| seq.insert_or_update_partially(*aabb, *leaf, *margin))
            .collect();

        let mut par_statuses = Vec::new();
        par.insert_or_update_batch_partially_parallel(&updates, &mut par_statuses);

        assert_eq!(seq_statuses, par_statuses);
        assert!(seq_statuses.contains(&BvhLeafUpdateStatus::Unchanged));
        assert!(seq_statuses.contains(&BvhLeafUpdateStatus::UpdatedInPlace));
        assert!(seq_statuses.contains(&BvhLeafUpdateStatus::Inserted));

        let mut w1 = BvhWorkspace::default();
        let mut w2 = BvhWorkspace::default();
        seq.refit(&mut w1);
        par.refit(&mut w2);
        for i in 0..1500u32 {
            assert_eq!(
                seq.leaf_node(i).map(|n| n.aabb()),
                par.leaf_node(i).map(|n| n.aabb()),
                "leaf {i} differs between sequential and parallel updates"
            );
        }
    }
}

/// The contracts rapier's `parallel`-off ≡ `parallel`-on guarantee rests on: when rapier
/// picks the parallel variant of one of these passes, it must get the exact result the
/// sequential one would have produced — otherwise a build with the feature and a build
/// without it are two different simulations.
#[cfg(feature = "parallel")]
#[cfg(all(feature = "dim3", feature = "f32"))]
mod parallel_matches_sequential {
    use crate::bounding_volume::Aabb;
    use crate::math::Vector;
    use crate::partitioning::{Bvh, BvhWorkspace};
    use alloc::vec::Vec;

    /// 17^3, past `refit_buffers_parallel`'s `SEQ_LEAF_THRESHOLD` so the parallel refit
    /// really splits instead of delegating to the sequential one.
    const SIDE: u32 = 17;
    const LEAVES: u32 = SIDE * SIDE * SIDE;

    /// A jittered grid whose cells overlap their neighbours: the tree needs real internal
    /// structure, and the traversal real pairs, for the comparisons to mean anything.
    fn scene() -> Bvh {
        let mut bvh = Bvh::new();
        for i in 0..LEAVES {
            let (x, y, z) = (i % SIDE, (i / SIDE) % SIDE, i / (SIDE * SIDE));
            let jitter = (i % 7) as f32 * 0.03;
            let mins = Vector::new(x as f32, y as f32, z as f32) * 1.0
                + Vector::new(jitter, jitter * 0.5, jitter * 0.25);
            bvh.insert(
                Aabb::new(mins.into(), (mins + Vector::new(1.5, 1.5, 1.5)).into()),
                i,
            );
        }
        bvh
    }

    /// Compares the node arrays themselves, not just the leaf AABBs: the refit rebuilds
    /// them in depth-first order, and that order decides the traversal order downstream.
    fn assert_same_nodes(seq: &Bvh, par: &Bvh) {
        assert_eq!(
            seq.nodes.len(),
            par.nodes.len(),
            "node array lengths differ after refit"
        );
        for (i, (a, b)) in seq.nodes.iter().zip(par.nodes.iter()).enumerate() {
            for (side, a, b) in [("left", &a.left, &b.left), ("right", &a.right, &b.right)] {
                assert_eq!(a.aabb(), b.aabb(), "node {i} ({side}) aabb differs");
                assert_eq!(
                    a.leaf_count(),
                    b.leaf_count(),
                    "node {i} ({side}) leaf count differs"
                );
                assert_eq!(
                    a.is_leaf(),
                    b.is_leaf(),
                    "node {i} ({side}) leafness differs"
                );
            }
        }
        for i in 0..LEAVES {
            assert_eq!(
                seq.leaf_node(i).map(|n| n.aabb()),
                par.leaf_node(i).map(|n| n.aabb()),
                "leaf {i} differs after refit"
            );
        }
    }

    /// [`Bvh::refit_parallel`] documents that "only the work distribution differs".
    #[test]
    fn refit_parallel_matches_sequential() {
        let (mut seq, mut par) = (scene(), scene());
        let (mut w1, mut w2) = (BvhWorkspace::default(), BvhWorkspace::default());
        seq.refit(&mut w1);
        par.refit_parallel(&mut w2);
        assert_same_nodes(&seq, &par);
    }

    /// Same, for the flag-preserving refit rapier runs in its deferred optimization pass.
    #[test]
    fn refit_without_resolve_parallel_matches_sequential() {
        let (mut seq, mut par) = (scene(), scene());
        let (mut w1, mut w2) = (BvhWorkspace::default(), BvhWorkspace::default());
        seq.optimize_incremental(&mut w1);
        par.optimize_incremental(&mut w2);
        seq.refit_without_resolve(&mut w1);
        par.refit_without_resolve_parallel(&mut w2);
        assert_same_nodes(&seq, &par);
    }

    /// [`Bvh::traverse_bvtt_single_tree_parallel`] documents that it returns the pairs "in
    /// the exact same order as the calls the sequential traversal would have made".
    #[test]
    fn bvtt_traversal_parallel_matches_sequential() {
        let mut bvh = scene();
        let mut workspace = BvhWorkspace::default();
        bvh.refit(&mut workspace);

        let mut sequential = Vec::new();
        bvh.traverse_bvtt_single_tree::<false>(&mut workspace, &mut |a, b| sequential.push((a, b)));
        let parallel = bvh.traverse_bvtt_single_tree_parallel::<false>();

        assert!(
            !sequential.is_empty(),
            "the scene must actually report pairs"
        );
        assert_eq!(
            sequential, parallel,
            "the parallel BVTT traversal reported a different pair sequence"
        );
    }
}
