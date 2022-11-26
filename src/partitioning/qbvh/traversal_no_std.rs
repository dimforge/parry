use crate::bounding_volume::{Aabb, SimdAabb};
use crate::math::Real;
use crate::partitioning::visitor::SimdSimultaneousVisitStatus;
use crate::partitioning::{
    GenericQbvh, QbvhStorage, SimdBestFirstVisitStatus, SimdBestFirstVisitor,
    SimdSimultaneousVisitor, SimdVisitStatus, SimdVisitor,
};
use crate::simd::SIMD_WIDTH;
use crate::utils::{Array1, WeightedValue};
use arrayvec::ArrayVec;
use num::Bounded;
use simba::simd::SimdBool;

use super::{IndexedData, NodeIndex, Qbvh};

impl<LeafData: IndexedData, Storage: QbvhStorage<LeafData>> GenericQbvh<LeafData, Storage> {
    /// Performs a depth-first traversal on the BVH.
    ///
    /// # Return
    ///
    /// Returns `false` if the traversal exited early, and `true` otherwise.
    pub fn traverse_depth_first(&self, visitor: &mut impl SimdVisitor<LeafData, SimdAabb>) -> bool {
        self.traverse_depth_first_node(visitor, 0)
    }

    /// Performs a depth-first traversal on the BVH starting at the given node.
    ///
    /// # Return
    ///
    /// Returns `false` if the traversal exited early, and `true` otherwise.
    pub fn traverse_depth_first_node(
        &self,
        visitor: &mut impl SimdVisitor<LeafData, SimdAabb>,
        curr_node: u32,
    ) -> bool {
        let node = &self.nodes[curr_node as usize];
        let leaf_data = if node.is_leaf() {
            Some(
                array![|ii| Some(&self.proxies.get_at(node.children[ii] as usize)?.data); SIMD_WIDTH],
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
                    if (bitmask & (1 << ii)) != 0 {
                        if !node.is_leaf() {
                            // Internal node, visit the child.
                            // Unfortunately, we have this check because invalid Aabbs
                            // return a hit as well.
                            if node.children[ii] as usize <= self.nodes.len() {
                                if !self.traverse_depth_first_node(visitor, node.children[ii]) {
                                    return false;
                                }
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
        self.traverse_best_first_node(visitor, 0, Real::MAX)
    }

    /// Performs a best-first-search on the BVH.
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

        let mut best_cost = init_cost;
        let mut result = None;

        // NOTE: a stack with 64 elements is enough for a depth-first search
        //       on a tree with up to about 4.000.000.000 triangles.
        //       See https://math.stackexchange.com/a/2739663 for the max
        //       stack depth on a depth-first search.
        let mut stack: ArrayVec<_, 64> = ArrayVec::new();
        stack.push(WeightedValue::new(start_node, -best_cost / 2.0));

        self.traverse_best_first_node_recursive(visitor, &mut stack, &mut best_cost, &mut result);
        result
    }

    fn traverse_best_first_node_recursive<BFS>(
        &self,
        visitor: &mut BFS,
        stack: &mut ArrayVec<WeightedValue<u32>, 64>,
        best_cost: &mut Real,
        best_result: &mut Option<(NodeIndex, BFS::Result)>,
    ) where
        BFS: SimdBestFirstVisitor<LeafData, SimdAabb>,
        BFS::Result: Clone, // Because we cannot move out of an array…
    {
        while let Some(entry) = stack.pop() {
            // NOTE: since we are not actually allowed to allocate, we can’t use the binary heap.
            //       So we are really just running a recursive depth-first traversal, with early
            //       exit on the current best cost.
            if -entry.cost >= *best_cost {
                continue;
            }

            let node = &self.nodes[entry.value as usize];
            let leaf_data = if node.is_leaf() {
                Some(
                    array![|ii| Some(&self.proxies.get_at(node.children[ii] as usize)?.data); SIMD_WIDTH],
                )
            } else {
                None
            };

            match visitor.visit(*best_cost, &node.simd_aabb, leaf_data) {
                SimdBestFirstVisitStatus::ExitEarly(result) => {
                    if result.is_some() {
                        *best_result = result.map(|r| (node.parent, r));
                        return;
                    }
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
                                if weights[ii] < *best_cost && results[ii].is_some() {
                                    // We found a leaf!
                                    if let Some(proxy) =
                                        self.proxies.get_at(node.children[ii] as usize)
                                    {
                                        *best_cost = weights[ii];
                                        *best_result =
                                            Some((proxy.node, results[ii].clone().unwrap()))
                                    }
                                }
                            } else {
                                // Internal node, visit the child.
                                // Unfortunately, we have this check because invalid Aabbs
                                // return a hit as well.
                                if (node.children[ii] as usize) < self.nodes.len() {
                                    let child_node =
                                        WeightedValue::new(node.children[ii], -weights[ii]);
                                    stack.push(child_node);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}
