use crate::bounding_volume::SimdAABB;
use crate::math::Real;
use crate::partitioning::{
    SimdBestFirstVisitStatus, SimdBestFirstVisitor, SimdVisitStatus, SimdVisitor,
};
use crate::simd::SIMD_WIDTH;
use crate::utils::WeightedValue;
use num::Bounded;
use simba::simd::SimdBool;
use std::collections::BinaryHeap;

use super::{IndexedData, NodeIndex, QBVH};

impl<T: IndexedData> QBVH<T> {
    /// Performs a depth-first traversal on the BVH.
    pub fn traverse_depth_first(&self, visitor: &mut impl SimdVisitor<T, SimdAABB>) {
        self.traverse_depth_first_with_stack(visitor, &mut Vec::new())
    }

    /// Performs a depth-first traversal on the BVH.
    pub fn traverse_depth_first_with_stack(
        &self,
        visitor: &mut impl SimdVisitor<T, SimdAABB>,
        stack: &mut Vec<u32>,
    ) {
        stack.clear();

        if !self.nodes.is_empty() {
            stack.push(0);
        }
        while let Some(entry) = stack.pop() {
            let node = self.nodes[entry as usize];
            let leaf_data = if node.leaf {
                Some(
                    array![|ii| Some(&self.proxies.get(node.children[ii] as usize)?.data); SIMD_WIDTH],
                )
            } else {
                None
            };

            match visitor.visit(&node.simd_aabb, leaf_data) {
                SimdVisitStatus::ExitEarly => {
                    return;
                }
                SimdVisitStatus::MaybeContinue(mask) => {
                    let bitmask = mask.bitmask();

                    for ii in 0..SIMD_WIDTH {
                        if (bitmask & (1 << ii)) != 0 {
                            if !node.leaf {
                                // Internal node, visit the child.
                                // Un fortunately, we have this check because invalid AABBs
                                // return a hit as well.
                                if node.children[ii] as usize <= self.nodes.len() {
                                    stack.push(node.children[ii]);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    /// Performs a best-first-search on the BVH.
    ///
    /// Returns the content of the leaf with the smallest associated cost, and a result of
    /// user-defined type.
    pub fn traverse_best_first<BFS>(&self, visitor: &mut BFS) -> Option<(NodeIndex, BFS::Result)>
    where
        BFS: SimdBestFirstVisitor<T, SimdAABB>,
        BFS::Result: Clone, // Because we cannot move out of an array…
    {
        if self.nodes.is_empty() {
            return None;
        }

        let mut queue: BinaryHeap<WeightedValue<u32>> = BinaryHeap::new();

        let mut best_cost = Real::max_value();
        let mut best_result = None;
        queue.push(WeightedValue::new(0, -best_cost / 2.0));

        while let Some(entry) = queue.pop() {
            if -entry.cost >= best_cost {
                // No BV left in the tree that has a lower cost than best_result
                break; // Solution found.
            }

            let node = self.nodes[entry.value as usize];
            let leaf_data = if node.leaf {
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
                            if node.leaf {
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
                                // Un fortunately, we have this check because invalid AABBs
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

    /// Retrieve all the data of the nodes with AABBs intersecting
    /// the given AABB:
    // FIXME: implement a visitor pattern to merge intersect_aabb
    // and intersect_ray into a single method.
    pub fn intersect_aabb(&self, aabb: &AABB, out: &mut Vec<T>) {
        if self.nodes.is_empty() {
            return;
        }

        // Special case for the root.
        let mut stack = vec![0u32];
        let simd_aabb = SimdAABB::splat(*aabb);
        while let Some(inode) = stack.pop() {
            let node = self.nodes[inode as usize];
            let intersections = node.simd_aabb.intersects(&simd_aabb);
            let bitmask = intersections.bitmask();

            for ii in 0..SIMD_WIDTH {
                if (bitmask & (1 << ii)) != 0 {
                    if node.leaf {
                        // We found a leaf!
                        // Unfortunately, invalid AABBs return a intersection as well.
                        if let Some(proxy) = self.proxies.get(node.children[ii] as usize) {
                            out.push(proxy.data);
                        }
                    } else {
                        // Internal node, visit the child.
                        // Unfortunately, we have this check because invalid AABBs
                        // return a intersection as well.
                        if node.children[ii] as usize <= self.nodes.len() {
                            stack.push(node.children[ii]);
                        }
                    }
                }
            }
        }
    }
}
