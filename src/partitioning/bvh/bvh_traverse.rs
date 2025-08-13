use super::BvhNode;
use crate::math::Real;
use crate::partitioning::Bvh;
use smallvec::SmallVec;

const TRAVERSAL_STACK_SIZE: usize = 32;

pub struct Leaves<'a, Check: Fn(&BvhNode) -> bool> {
    tree: &'a Bvh,
    next: Option<&'a BvhNode>,
    stack: SmallVec<[&'a BvhNode; TRAVERSAL_STACK_SIZE]>,
    check: Check,
}

impl<'a, Check: Fn(&BvhNode) -> bool> Iterator for Leaves<'a, Check> {
    type Item = u32;
    fn next(&mut self) -> Option<Self::Item> {
        loop {
            if self.next.is_none() {
                self.next = self.stack.pop();
            }

            let node = self.next.take()?;

            if node.is_leaf() {
                if (self.check)(node) {
                    return Some(node.children);
                } else {
                    continue;
                }
            }
            
            let children = &self.tree.nodes[node.children as usize];
            let left = &children.left;
            let right = &children.right;

            if (self.check)(left) {
                self.next = Some(left);
            }

            if (self.check)(right) {
                if self.next.is_none() {
                    self.next = Some(right);
                } else {
                    self.stack.push(right);
                }
            }
        }
    }
}

/// Cost associated to a BVH leaf during best-first traversal.
pub trait BvhLeafCost {
    /// The cost value associated to the leaf.
    ///
    /// Best-first searches for the leaf with the lowest cost.
    fn cost(&self) -> Real;
}

impl BvhLeafCost for Real {
    #[inline(always)]
    fn cost(&self) -> Real {
        *self
    }
}

impl<T> BvhLeafCost for (Real, T) {
    #[inline(always)]
    fn cost(&self) -> Real {
        self.0
    }
}

impl Bvh {
    /// Iterates through the leaves, in depth-first order.
    ///
    /// The `check_node` closure is called on every traversed node. If it returns `false` then the
    /// node and all its descendants won’t be iterated on. This is useful for pruning whole
    /// sub-trees based on a geometric predicate on the node’s AABB.
    ///
    /// See also the [`Bvh::traverse`] function which is slightly less convenient since it doesn’t
    /// rely on the iterator system, but takes a closure that implements [`FnMut`] instead of [`Fn`].
    pub fn leaves<F: Fn(&BvhNode) -> bool>(&self, check_node: F) -> Leaves<'_, F> {
        if let Some(root) = self.nodes.first() {
            let mut stack = SmallVec::default();
            if root.right.leaf_count() > 0 {
                stack.push(&root.right);
            }

            Leaves {
                tree: self,
                next: Some(&root.left),
                stack,
                check: check_node,
            }
        } else {
            Leaves {
                tree: self,
                next: None,
                stack: Default::default(),
                check: check_node,
            }
        }
    }
}

/// Controls the execution flowo of [`Bvh::traverse`].
pub enum TraversalAction {
    /// The traversal will continue on the children of the tested node.
    Continue,
    /// The traversal will skip all descendants of the tested node.
    Prune,
    /// The traversal will exit immediately.
    EarlyExit,
}

impl Bvh {
    #[inline(always)]
    pub(crate) fn traversal_stack() -> SmallVec<[u32; 32]> {
        Default::default()
    }

    /// Traverse the tree in depth-first order.
    ///
    /// The `check_node` closure is called on every traversed node. The returned [`TraversalAction`]
    /// controls whether a given node (and all its subtree) needs to be traversed, skipped, or if
    /// the traversal needs to exit immediately (for example if you were looking for only one
    /// particular node).
    ///
    /// See also the [`Bvh::leaves`] iterator which is a more convenient way of traversing the tree,
    /// but is slightly limited in terms of node checking. In particular the closure
    /// given to [`Bvh::traverse`] is mutable can hold any state so its check can depend on previous
    /// execution of itself during the traversal.
    pub fn traverse(&self, mut check_node: impl FnMut(&BvhNode) -> TraversalAction) {
        let mut stack = Self::traversal_stack();
        let mut curr_id = 0;

        if self.nodes.is_empty() {
            return;
        } else if self.nodes[0].right.leaf_count() == 0 {
            // Special case for partial root.
            let _ = check_node(&self.nodes[0].left);
            return;
        }

        loop {
            let node = &self.nodes[curr_id as usize];
            let left = &node.left;
            let right = &node.right;
            let go_left = match check_node(left) {
                TraversalAction::Continue => !left.is_leaf(),
                TraversalAction::Prune => false,
                TraversalAction::EarlyExit => return,
            };
            let go_right = match check_node(right) {
                TraversalAction::Continue => !right.is_leaf(),
                TraversalAction::Prune => false,
                TraversalAction::EarlyExit => return,
            };

            match (go_left, go_right) {
                (true, true) => {
                    curr_id = left.children;
                    stack.push(right.children);
                }
                (true, false) => curr_id = left.children,
                (false, true) => curr_id = right.children,
                (false, false) => {
                    let Some(next) = stack.pop() else {
                        return;
                    };
                    curr_id = next;
                }
            }
        }
    }

    /// Find the leaf that minimizes their associated cost.
    pub fn find_best<L: BvhLeafCost>(
        &self,
        max_cost: Real,
        aabb_cost: impl Fn(&BvhNode, Real) -> Real,
        leaf_cost: impl Fn(u32, Real) -> Option<L>,
    ) -> Option<(u32, L)> {
        // A stack with 32 elements should be more than enough in most cases.
        let mut stack = Self::traversal_stack();
        let mut best_val = None;
        let mut best_cost = max_cost;
        let mut best_id = u32::MAX;
        let mut curr_id = 0;

        if self.nodes.is_empty() {
            return None;
        } else if self.nodes[0].right.leaf_count() == 0 {
            // Special case for partial root.
            let leaf = &self.nodes[0].left;
            if aabb_cost(leaf, max_cost) < max_cost {
                let cost = leaf_cost(leaf.children, best_cost)?;
                return (cost.cost() < max_cost).then_some((leaf.children, cost));
            } else {
                return None;
            }
        }

        loop {
            let node = &self.nodes[curr_id as usize];
            let mut left = &node.left;
            let mut right = &node.right;

            let mut left_score = aabb_cost(left, best_cost);
            let mut right_score = aabb_cost(right, best_cost);

            if left_score > right_score {
                core::mem::swap(&mut left_score, &mut right_score);
                core::mem::swap(&mut left, &mut right);
            }

            let mut found_next = false;
            if left_score < best_cost && left_score != Real::MAX {
                if left.is_leaf() {
                    if let Some(primitive_val) = leaf_cost(left.children, best_cost) {
                        let primitive_score = primitive_val.cost();
                        if primitive_score < best_cost {
                            best_val = Some(primitive_val);
                            best_cost = primitive_score;
                            best_id = left.children;
                        }
                    }
                } else {
                    curr_id = left.children;
                    found_next = true;
                }
            }

            if right_score < best_cost && right_score != Real::MAX {
                if right.is_leaf() {
                    if let Some(primitive_val) = leaf_cost(right.children, best_cost) {
                        let primitive_score = primitive_val.cost();
                        if primitive_score < best_cost {
                            best_val = Some(primitive_val);
                            best_cost = primitive_score;
                            best_id = right.children;
                        }
                    }
                } else if found_next {
                    stack.push(right.children);
                } else {
                    curr_id = right.children;
                    found_next = true;
                }
            }

            if !found_next {
                if let Some(next) = stack.pop() {
                    curr_id = next;
                } else {
                    return best_val.map(|val| (best_id, val));
                }
            }
        }
    }
}
