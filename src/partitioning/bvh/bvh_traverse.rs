use super::BvhNode;
use crate::math::Real;
use crate::partitioning::Bvh;
use crate::query::{RayIntersection, ShapeCastHit};
use smallvec::SmallVec;

const TRAVERSAL_STACK_SIZE: usize = 32;

pub trait NodeCheck {
    fn check(&self, node: &BvhNode) -> bool;
}

impl<F: Fn(&BvhNode) -> bool> NodeCheck for F {
    #[inline(always)]
    fn check(&self, node: &BvhNode) -> bool {
        self(node)
    }
}

impl NodeCheck for () {
    #[inline(always)]
    fn check(&self, _: &BvhNode) -> bool {
        true
    }
}

pub struct Leaves<'a, Check: NodeCheck> {
    tree: &'a Bvh,
    next: Option<&'a BvhNode>,
    stack: SmallVec<[&'a BvhNode; TRAVERSAL_STACK_SIZE]>,
    check: Check,
}

impl<'a, Check: NodeCheck> Iterator for Leaves<'a, Check> {
    type Item = u32;
    fn next(&mut self) -> Option<Self::Item> {
        loop {
            if self.next.is_none() {
                self.next = self.stack.pop();
            }

            let node = self.next.take()?;

            if node.is_leaf() {
                return Some(node.children);
            }

            let children = &self.tree.nodes[node.children as usize];
            let left = &children.left;
            let right = &children.right;

            if self.check.check(left) {
                self.next = Some(left);
            }

            if self.check.check(right) {
                if self.next.is_none() {
                    self.next = Some(right);
                } else {
                    self.stack.push(right);
                }
            }
        }
    }
}

pub trait AabbCost {
    fn cost(&self, node: &BvhNode, best_so_far: Real) -> Real;
}

pub trait LeafCost {
    type CostValue: LeafCostValue;
    fn cost(&self, leaf: u32, best_so_far: Real) -> Option<Self::CostValue>;
}

impl<F: Fn(&BvhNode, Real) -> Real> AabbCost for F {
    #[inline(always)]
    fn cost(&self, node: &BvhNode, best_so_far: Real) -> Real {
        self(node, best_so_far)
    }
}

impl<C: LeafCostValue, F: Fn(u32, Real) -> Option<C>> LeafCost for F {
    type CostValue = C;

    #[inline(always)]
    fn cost(&self, leaf: u32, best_so_far: Real) -> Option<Self::CostValue> {
        self(leaf, best_so_far)
    }
}

pub trait LeafCostValue {
    fn cost(&self) -> Real;
}

impl LeafCostValue for Real {
    #[inline(always)]
    fn cost(&self) -> Real {
        *self
    }
}

impl<T> LeafCostValue for (Real, T) {
    #[inline(always)]
    fn cost(&self) -> Real {
        self.0
    }
}

// TODO: move this to ray.rs
impl LeafCostValue for RayIntersection {
    #[inline]
    fn cost(&self) -> Real {
        self.time_of_impact
    }
}

impl LeafCostValue for ShapeCastHit {
    #[inline]
    fn cost(&self) -> Real {
        self.time_of_impact
    }
}

impl Bvh {
    pub fn leaves<F: NodeCheck>(&self, check: F) -> Leaves<F> {
        if let Some(root) = self.nodes.first() {
            let mut stack = SmallVec::default();
            if root.right.leaf_count() > 0 {
                stack.push(&root.right);
            }

            Leaves {
                tree: self,
                next: Some(&root.left),
                stack,
                check,
            }
        } else {
            Leaves {
                tree: self,
                next: None,
                stack: Default::default(),
                check,
            }
        }
    }
}

pub enum TraversalAction {
    Continue,
    Prune,
    EarlyExit,
}

pub trait NodeVisitor {
    fn visit(&mut self, node: &BvhNode) -> TraversalAction;
}

impl<F: FnMut(&BvhNode) -> TraversalAction> NodeVisitor for F {
    #[inline(always)]
    fn visit(&mut self, node: &BvhNode) -> TraversalAction {
        self(node)
    }
}

impl Bvh {
    #[inline(always)]
    pub(crate) fn traversal_stack() -> SmallVec<[u32; 32]> {
        Default::default()
    }

    pub fn traverse(&self, mut check_node: impl NodeVisitor) {
        let mut stack = Self::traversal_stack();
        let mut curr_id = 0;

        if self.nodes.is_empty() {
            return;
        } else if self.nodes[0].right.leaf_count() == 0 {
            // Special case for partial root.
            let _ = check_node.visit(&self.nodes[0].left);
            return;
        }

        loop {
            let node = &self.nodes[curr_id as usize];
            let left = &node.left;
            let right = &node.right;
            let go_left = match check_node.visit(left) {
                TraversalAction::Continue => !left.is_leaf(),
                TraversalAction::Prune => false,
                TraversalAction::EarlyExit => return,
            };
            let go_right = match check_node.visit(right) {
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
    pub fn find_best<L: LeafCost>(
        &self,
        max_cost: Real,
        aabb_cost: impl AabbCost, // impl Fn(&BvhNode, Real) -> Real,
        leaf_cost: L,             // impl Fn(u32, Real) -> L,
    ) -> Option<(u32, L::CostValue)> {
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
            if aabb_cost.cost(leaf, max_cost) < max_cost {
                let cost = leaf_cost.cost(leaf.children, best_cost)?;
                return (cost.cost() < max_cost).then_some((leaf.children, cost));
            } else {
                return None;
            }
        }

        loop {
            let node = &self.nodes[curr_id as usize];
            let mut left = &node.left;
            let mut right = &node.right;

            let mut left_score = aabb_cost.cost(left, best_cost);
            let mut right_score = aabb_cost.cost(right, best_cost);

            if left_score > right_score {
                core::mem::swap(&mut left_score, &mut right_score);
                core::mem::swap(&mut left, &mut right);
            }

            let mut found_next = false;
            if left_score < best_cost && left_score != Real::MAX {
                if left.is_leaf() {
                    if let Some(primitive_val) = leaf_cost.cost(left.children, best_cost) {
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
                    if let Some(primitive_val) = leaf_cost.cost(right.children, best_cost) {
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
