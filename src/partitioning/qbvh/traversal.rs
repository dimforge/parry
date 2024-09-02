#![allow(clippy::needless_range_loop)] // This tends to make the traversal code much more verbose that necessary.

use crate::bounding_volume::{Aabb, SimdAabb};
use crate::math::Real;
use crate::partitioning::visitor::{SimdSimultaneousVisitStatus, SimdVisitorWithContext};
use crate::partitioning::{
    Qbvh, SimdBestFirstVisitStatus, SimdBestFirstVisitor, SimdSimultaneousVisitor, SimdVisitStatus,
    SimdVisitor,
};
use crate::simd::SIMD_WIDTH;
use crate::utils::WeightedValue;
use num::Bounded;
use simba::simd::SimdBool;
use std::collections::BinaryHeap;
#[cfg(feature = "parallel")]
use {
    crate::partitioning::{ParallelSimdSimultaneousVisitor, ParallelSimdVisitor},
    arrayvec::ArrayVec,
    rayon::prelude::*,
    std::sync::atomic::{AtomicBool, Ordering as AtomicOrdering},
};

use super::{IndexedData, NodeIndex};

impl<LeafData: IndexedData> Qbvh<LeafData> {
    /// Performs a depth-first traversal on the BVH.
    ///
    /// # Return
    ///
    /// Returns `false` if the traversal exited early, and `true` otherwise.
    pub fn traverse_depth_first(&self, visitor: &mut impl SimdVisitor<LeafData, SimdAabb>) -> bool {
        self.traverse_depth_first_node(visitor, 0)
    }

    /// Performs a depth-first traversal on the BVH, starting at the given node.
    ///
    /// # Return
    ///
    /// Returns `false` if the traversal exited early, and `true` otherwise.
    pub fn traverse_depth_first_node(
        &self,
        visitor: &mut impl SimdVisitor<LeafData, SimdAabb>,
        start_node: u32,
    ) -> bool {
        self.traverse_depth_first_node_with_stack(visitor, &mut Vec::new(), start_node)
    }

    /// Performs a depth-first traversal on the BVH.
    ///
    /// # Return
    ///
    /// Returns `false` if the traversal exited early, and `true` otherwise.
    pub fn traverse_depth_first_with_stack(
        &self,
        visitor: &mut impl SimdVisitor<LeafData, SimdAabb>,
        stack: &mut Vec<u32>,
    ) -> bool {
        self.traverse_depth_first_node_with_stack(visitor, stack, 0)
    }

    /// Performs a depth-first traversal on the BVH.
    ///
    /// # Return
    ///
    /// Returns `false` if the traversal exited early, and `true` otherwise.
    pub fn traverse_depth_first_node_with_stack(
        &self,
        visitor: &mut impl SimdVisitor<LeafData, SimdAabb>,
        stack: &mut Vec<u32>,
        start_node: u32,
    ) -> bool {
        stack.clear();

        if !self.nodes.is_empty() {
            stack.push(start_node);
        }
        while let Some(entry) = stack.pop() {
            let node = &self.nodes[entry as usize];
            let leaf_data = if node.is_leaf() {
                Some(
                    array![|ii| Some(&self.proxies.get(node.children[ii] as usize)?.data); SIMD_WIDTH],
                )
            } else {
                None
            };

            match visitor.visit(&node.simd_aabb, leaf_data) {
                SimdVisitStatus::ExitEarly => {
                    return false;
                }
                SimdVisitStatus::MaybeContinue(mask) => {
                    let bitmask = mask.bitmask();

                    for ii in 0..SIMD_WIDTH {
                        if (bitmask & (1 << ii)) != 0 && !node.is_leaf() {
                            // Internal node, visit the child.
                            // Un fortunately, we have this check because invalid Aabbs
                            // return a hit as well.
                            if node.children[ii] as usize <= self.nodes.len() {
                                stack.push(node.children[ii]);
                            }
                        }
                    }
                }
            }
        }

        true
    }

    /// Performs a depth-first traversal on the BVH. Passes a context from the
    /// parent to the children.
    ///
    /// # Return
    ///
    /// Returns `false` if the traversal exited early, and `true` otherwise.
    pub fn traverse_depth_first_with_context<Context: Clone>(
        &self,
        visitor: &mut impl SimdVisitorWithContext<LeafData, SimdAabb, Context>,
        context: Context,
    ) -> bool {
        self.traverse_depth_first_node_with_stack_and_context(visitor, &mut Vec::new(), 0, context)
    }

    /// Performs a depth-first traversal on the BVH and propagates a context down,
    /// from the root to each of its descendants. The context can be modified
    /// during the query.
    ///
    /// # Return
    ///
    /// Returns `false` if the traversal exited early, and `true` otherwise.
    pub fn traverse_depth_first_node_with_stack_and_context<Context: Clone>(
        &self,
        visitor: &mut impl SimdVisitorWithContext<LeafData, SimdAabb, Context>,
        stack: &mut Vec<(u32, Context)>,
        start_node: u32,
        context: Context,
    ) -> bool {
        stack.clear();

        if !self.nodes.is_empty() {
            stack.push((start_node, context));
        }
        while let Some((entry, context)) = stack.pop() {
            let node = &self.nodes[entry as usize];
            let leaf_data = if node.is_leaf() {
                Some(
                    array![|ii| Some(&self.proxies.get(node.children[ii] as usize)?.data); SIMD_WIDTH],
                )
            } else {
                None
            };

            let (visit_result, contexts) = visitor.visit(&node.simd_aabb, leaf_data, context);
            match visit_result {
                SimdVisitStatus::ExitEarly => {
                    return false;
                }
                SimdVisitStatus::MaybeContinue(mask) => {
                    let bitmask = mask.bitmask();

                    for ii in 0..SIMD_WIDTH {
                        if (bitmask & (1 << ii)) != 0 && !node.is_leaf() {
                            // Internal node, visit the child.
                            // Un fortunately, we have this check because invalid Aabbs
                            // return a hit as well.
                            if node.children[ii] as usize <= self.nodes.len() {
                                stack.push((node.children[ii], contexts[ii].clone()));
                            }
                        }
                    }
                }
            }
        }

        true
    }

    /// Performs a best-first-search on the BVH.
    ///
    /// Returns the content of the leaf with the smallest associated cost, and a result of
    /// user-defined type.
    pub fn traverse_best_first<BFS>(&self, visitor: &mut BFS) -> Option<(NodeIndex, BFS::Result)>
    where
        BFS: SimdBestFirstVisitor<LeafData, SimdAabb>,
        BFS::Result: Clone, // Because we cannot move out of an array…
    {
        self.traverse_best_first_node(visitor, 0, Real::max_value())
    }

    /// Performs a best-first-search on the BVH, starting at the given node.
    ///
    /// Returns the content of the leaf with the smallest associated cost, and a result of
    /// user-defined type.
    pub fn traverse_best_first_node<BFS>(
        &self,
        visitor: &mut BFS,
        start_node: u32,
        init_cost: Real,
    ) -> Option<(NodeIndex, BFS::Result)>
    where
        BFS: SimdBestFirstVisitor<LeafData, SimdAabb>,
        BFS::Result: Clone, // Because we cannot move out of an array…
    {
        if self.nodes.is_empty() {
            return None;
        }

        let mut queue: BinaryHeap<WeightedValue<u32>> = BinaryHeap::new();

        let mut best_cost = init_cost;
        let mut best_result = None;
        queue.push(WeightedValue::new(start_node, -best_cost / 2.0));

        while let Some(entry) = queue.pop() {
            if -entry.cost >= best_cost {
                // No BV left in the tree that has a lower cost than best_result
                break; // Solution found.
            }

            let node = &self.nodes[entry.value as usize];
            let leaf_data = if node.is_leaf() {
                Some(
                    array![|ii| Some(&self.proxies.get(node.children[ii] as usize)?.data); SIMD_WIDTH],
                )
            } else {
                None
            };

            match visitor.visit(best_cost, &node.simd_aabb, leaf_data) {
                SimdBestFirstVisitStatus::ExitEarly(result) => {
                    return result.map(|r| (node.parent, r)).or(best_result);
                }
                SimdBestFirstVisitStatus::MaybeContinue {
                    weights,
                    mask,
                    results,
                } => {
                    let bitmask = mask.bitmask();
                    let weights: [Real; SIMD_WIDTH] = weights.into();

                    for ii in 0..SIMD_WIDTH {
                        if (bitmask & (1 << ii)) != 0 {
                            if node.is_leaf() {
                                if weights[ii] < best_cost && results[ii].is_some() {
                                    // We found a leaf!
                                    if let Some(proxy) =
                                        self.proxies.get(node.children[ii] as usize)
                                    {
                                        best_cost = weights[ii];
                                        best_result =
                                            Some((proxy.node, results[ii].clone().unwrap()))
                                    }
                                }
                            } else {
                                // Internal node, visit the child.
                                // Un fortunately, we have this check because invalid Aabbs
                                // return a hit as well.
                                if (node.children[ii] as usize) < self.nodes.len() {
                                    queue.push(WeightedValue::new(node.children[ii], -weights[ii]));
                                }
                            }
                        }
                    }
                }
            }
        }

        best_result
    }

    /// Retrieve all the data of the nodes with Aabbs intersecting
    /// the given Aabb:
    // TODO: implement a visitor pattern to merge intersect_aabb
    // and intersect_ray into a single method.
    pub fn intersect_aabb(&self, aabb: &Aabb, out: &mut Vec<LeafData>) {
        if self.nodes.is_empty() {
            return;
        }

        // Special case for the root.
        let mut stack = vec![0u32];
        let simd_aabb = SimdAabb::splat(*aabb);
        while let Some(inode) = stack.pop() {
            let node = &self.nodes[inode as usize];
            let intersections = node.simd_aabb.intersects(&simd_aabb);
            let bitmask = intersections.bitmask();

            for ii in 0..SIMD_WIDTH {
                if (bitmask & (1 << ii)) != 0 {
                    if node.is_leaf() {
                        // We found a leaf!
                        // Unfortunately, invalid Aabbs return a intersection as well.
                        if let Some(proxy) = self.proxies.get(node.children[ii] as usize) {
                            out.push(proxy.data);
                        }
                    } else {
                        // Internal node, visit the child.
                        // Unfortunately, we have this check because invalid Aabbs
                        // return a intersection as well.
                        if node.children[ii] as usize <= self.nodes.len() {
                            stack.push(node.children[ii]);
                        }
                    }
                }
            }
        }
    }

    /// Performs a simultaneous traversal of two Qbvh.
    pub fn traverse_bvtt<LeafData2: IndexedData>(
        &self,
        qbvh2: &Qbvh<LeafData2>,
        visitor: &mut impl SimdSimultaneousVisitor<LeafData, LeafData2, SimdAabb>,
    ) {
        self.traverse_bvtt_with_stack(qbvh2, visitor, &mut Vec::new())
    }

    /// Performs a simultaneous traversal of two Qbvh.
    pub fn traverse_bvtt_with_stack<LeafData2: IndexedData>(
        &self,
        qbvh2: &Qbvh<LeafData2>,
        visitor: &mut impl SimdSimultaneousVisitor<LeafData, LeafData2, SimdAabb>,
        stack: &mut Vec<(u32, u32)>,
    ) {
        let qbvh1 = self;
        stack.clear();

        if !qbvh1.nodes.is_empty() && !qbvh2.nodes.is_empty() {
            stack.push((0, 0));
        }

        while let Some(entry) = stack.pop() {
            let node1 = &qbvh1.nodes[entry.0 as usize];
            let node2 = &qbvh2.nodes[entry.1 as usize];

            let leaf_data1 = if node1.is_leaf() {
                Some(
                    array![|ii| Some(&qbvh1.proxies.get(node1.children[ii] as usize)?.data); SIMD_WIDTH],
                )
            } else {
                None
            };

            let leaf_data2 = if node2.is_leaf() {
                Some(
                    array![|ii| Some(&qbvh2.proxies.get(node2.children[ii] as usize)?.data); SIMD_WIDTH],
                )
            } else {
                None
            };

            match visitor.visit(&node1.simd_aabb, leaf_data1, &node2.simd_aabb, leaf_data2) {
                SimdSimultaneousVisitStatus::ExitEarly => {
                    return;
                }
                SimdSimultaneousVisitStatus::MaybeContinue(mask) => {
                    match (node1.is_leaf(), node2.is_leaf()) {
                        (true, true) => { /* Can’t go deeper. */ }
                        (true, false) => {
                            let mut bitmask = 0;
                            for ii in 0..SIMD_WIDTH {
                                bitmask |= mask[ii].bitmask();
                            }

                            for jj in 0..SIMD_WIDTH {
                                if (bitmask & (1 << jj)) != 0
                                    && node2.children[jj] as usize <= qbvh2.nodes.len()
                                {
                                    stack.push((entry.0, node2.children[jj]));
                                }
                            }
                        }
                        (false, true) => {
                            for ii in 0..SIMD_WIDTH {
                                let bitmask = mask[ii].bitmask();

                                if bitmask != 0 && node1.children[ii] as usize <= qbvh1.nodes.len()
                                {
                                    stack.push((node1.children[ii], entry.1));
                                }
                            }
                        }
                        (false, false) => {
                            for ii in 0..SIMD_WIDTH {
                                let bitmask = mask[ii].bitmask();

                                for jj in 0..SIMD_WIDTH {
                                    if (bitmask & (1 << jj)) != 0
                                        && node1.children[ii] as usize <= qbvh1.nodes.len()
                                        && node2.children[jj] as usize <= qbvh2.nodes.len()
                                    {
                                        stack.push((node1.children[ii], node2.children[jj]));
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    /// Performs a simultaneous traversal of two Qbvh.
    pub fn traverse_modified_bvtt<LeafData2: IndexedData>(
        &self,
        qbvh2: &Qbvh<LeafData2>,
        visitor: &mut impl SimdSimultaneousVisitor<LeafData, LeafData2, SimdAabb>,
    ) {
        self.traverse_modified_bvtt_with_stack(qbvh2, visitor, &mut Vec::new())
    }

    /// Performs a simultaneous traversal of two Qbvh.
    pub fn traverse_modified_bvtt_with_stack<LeafData2: IndexedData>(
        &self,
        qbvh2: &Qbvh<LeafData2>,
        visitor: &mut impl SimdSimultaneousVisitor<LeafData, LeafData2, SimdAabb>,
        stack: &mut Vec<(u32, u32)>,
    ) {
        let qbvh1 = self;
        stack.clear();

        if !qbvh1.nodes.is_empty() && !qbvh2.nodes.is_empty() && qbvh1.nodes[0].is_changed() {
            stack.push((0, 0));
        }

        while let Some(entry) = stack.pop() {
            let node1 = &qbvh1.nodes[entry.0 as usize];
            let node2 = &qbvh2.nodes[entry.1 as usize];

            if !node1.is_changed() {
                continue;
            }

            let leaf_data1 = if node1.is_leaf() {
                Some(
                    array![|ii| Some(&qbvh1.proxies.get(node1.children[ii] as usize)?.data); SIMD_WIDTH],
                )
            } else {
                None
            };

            let leaf_data2 = if node2.is_leaf() {
                Some(
                    array![|ii| Some(&qbvh2.proxies.get(node2.children[ii] as usize)?.data); SIMD_WIDTH],
                )
            } else {
                None
            };

            match visitor.visit(&node1.simd_aabb, leaf_data1, &node2.simd_aabb, leaf_data2) {
                SimdSimultaneousVisitStatus::ExitEarly => {
                    return;
                }
                SimdSimultaneousVisitStatus::MaybeContinue(mask) => {
                    match (node1.is_leaf(), node2.is_leaf()) {
                        (true, true) => { /* Can’t go deeper. */ }
                        (true, false) => {
                            let mut bitmask = 0;
                            for ii in 0..SIMD_WIDTH {
                                bitmask |= mask[ii].bitmask();
                            }

                            for jj in 0..SIMD_WIDTH {
                                if (bitmask & (1 << jj)) != 0
                                    && node2.children[jj] as usize <= qbvh2.nodes.len()
                                {
                                    stack.push((entry.0, node2.children[jj]));
                                }
                            }
                        }
                        (false, true) => {
                            for ii in 0..SIMD_WIDTH {
                                let bitmask = mask[ii].bitmask();

                                if bitmask != 0 && node1.children[ii] as usize <= qbvh1.nodes.len()
                                {
                                    stack.push((node1.children[ii], entry.1));
                                }
                            }
                        }
                        (false, false) => {
                            for ii in 0..SIMD_WIDTH {
                                let bitmask = mask[ii].bitmask();

                                for jj in 0..SIMD_WIDTH {
                                    if (bitmask & (1 << jj)) != 0
                                        && node1.children[ii] as usize <= qbvh1.nodes.len()
                                        && node2.children[jj] as usize <= qbvh2.nodes.len()
                                    {
                                        stack.push((node1.children[ii], node2.children[jj]));
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

#[cfg(feature = "parallel")]
impl<LeafData: IndexedData + Sync> Qbvh<LeafData> {
    /// Performs a depth-first traversal of two Qbvh using
    /// parallelism internally for better performances with large tree.
    pub fn traverse_depth_first_parallel(&self, visitor: &impl ParallelSimdVisitor<LeafData>) {
        if !self.nodes.is_empty() {
            let exit_early = AtomicBool::new(false);
            self.traverse_depth_first_node_parallel(visitor, &exit_early, 0);
        }
    }

    /// Runs a parallel depth-first traversal of the sub-tree starting at the given node.
    pub fn traverse_depth_first_node_parallel(
        &self,
        visitor: &impl ParallelSimdVisitor<LeafData>,
        exit_early: &AtomicBool,
        entry: u32,
    ) {
        if exit_early.load(AtomicOrdering::Relaxed) {
            return;
        }

        let mut stack: ArrayVec<u32, SIMD_WIDTH> = ArrayVec::new();
        let node = &self.nodes[entry as usize];
        let leaf_data = if node.is_leaf() {
            Some(array![|ii| Some(&self.proxies.get(node.children[ii] as usize)?.data); SIMD_WIDTH])
        } else {
            None
        };

        match visitor.visit(entry, node, leaf_data) {
            SimdVisitStatus::ExitEarly => {
                exit_early.store(true, AtomicOrdering::Relaxed);
                return;
            }
            SimdVisitStatus::MaybeContinue(mask) => {
                let bitmask = mask.bitmask();

                for ii in 0..SIMD_WIDTH {
                    if (bitmask & (1 << ii)) != 0 {
                        if !node.is_leaf() {
                            // Internal node, visit the child.
                            // Un fortunately, we have this check because invalid Aabbs
                            // return a hit as well.
                            if node.children[ii] as usize <= self.nodes.len() {
                                stack.push(node.children[ii]);
                            }
                        }
                    }
                }
            }
        }

        stack
            .as_slice()
            .par_iter()
            .copied()
            .for_each(|entry| self.traverse_depth_first_node_parallel(visitor, exit_early, entry));
    }

    /// Performs a simultaneous traversal of two Qbvh using
    /// parallelism internally for better performances with large tree.
    pub fn traverse_bvtt_parallel<
        LeafData2: IndexedData + Sync,
        Visitor: ParallelSimdSimultaneousVisitor<LeafData, LeafData2>,
    >(
        &self,
        qbvh2: &Qbvh<LeafData2>,
        visitor: &Visitor,
    ) {
        if !self.nodes.is_empty() && !qbvh2.nodes.is_empty() {
            let exit_early = AtomicBool::new(false);
            self.traverse_bvtt_node_parallel(
                qbvh2,
                visitor,
                &exit_early,
                Visitor::Data::default(),
                (0, 0),
            );
        }
    }

    /// Runs a parallel simultaneous traversal of the sub-tree starting at the given nodes.
    pub fn traverse_bvtt_node_parallel<
        LeafData2: IndexedData + Sync,
        Visitor: ParallelSimdSimultaneousVisitor<LeafData, LeafData2>,
    >(
        &self,
        qbvh2: &Qbvh<LeafData2>,
        visitor: &Visitor,
        exit_early: &AtomicBool,
        data: Visitor::Data,
        entry: (u32, u32),
    ) {
        if exit_early.load(AtomicOrdering::Relaxed) {
            return;
        }

        let qbvh1 = self;
        let node1 = &qbvh1.nodes[entry.0 as usize];
        let node2 = &qbvh2.nodes[entry.1 as usize];

        const SQUARE_SIMD_WIDTH: usize = SIMD_WIDTH * SIMD_WIDTH;
        let mut stack: ArrayVec<(u32, u32), SQUARE_SIMD_WIDTH> = ArrayVec::new();

        let leaf_data1 = if node1.is_leaf() {
            Some(
                array![|ii| Some(&qbvh1.proxies.get(node1.children[ii] as usize)?.data); SIMD_WIDTH],
            )
        } else {
            None
        };

        let leaf_data2 = if node2.is_leaf() {
            Some(
                array![|ii| Some(&qbvh2.proxies.get(node2.children[ii] as usize)?.data); SIMD_WIDTH],
            )
        } else {
            None
        };

        let (status, data) = visitor.visit(
            entry.0, &node1, leaf_data1, entry.1, &node2, leaf_data2, data,
        );

        match status {
            SimdSimultaneousVisitStatus::ExitEarly => {
                exit_early.store(true, AtomicOrdering::Relaxed);
                return;
            }
            SimdSimultaneousVisitStatus::MaybeContinue(mask) => {
                match (node1.is_leaf(), node2.is_leaf()) {
                    (true, true) => { /* Can’t go deeper. */ }
                    (true, false) => {
                        let mut bitmask = 0;
                        for ii in 0..SIMD_WIDTH {
                            bitmask |= mask[ii].bitmask();
                        }

                        for jj in 0..SIMD_WIDTH {
                            if (bitmask & (1 << jj)) != 0 {
                                if node2.children[jj] as usize <= qbvh2.nodes.len() {
                                    stack.push((entry.0, node2.children[jj]));
                                }
                            }
                        }
                    }
                    (false, true) => {
                        for ii in 0..SIMD_WIDTH {
                            let bitmask = mask[ii].bitmask();

                            if bitmask != 0 {
                                if node1.children[ii] as usize <= qbvh1.nodes.len() {
                                    stack.push((node1.children[ii], entry.1));
                                }
                            }
                        }
                    }
                    (false, false) => {
                        for ii in 0..SIMD_WIDTH {
                            let bitmask = mask[ii].bitmask();

                            for jj in 0..SIMD_WIDTH {
                                if (bitmask & (1 << jj)) != 0 {
                                    if node1.children[ii] as usize <= qbvh1.nodes.len()
                                        && node2.children[jj] as usize <= qbvh2.nodes.len()
                                    {
                                        stack.push((node1.children[ii], node2.children[jj]));
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        stack.as_slice().par_iter().copied().for_each(|entry| {
            self.traverse_bvtt_node_parallel(qbvh2, visitor, exit_early, data, entry)
        });
    }
}
