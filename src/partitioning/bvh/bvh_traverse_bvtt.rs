use super::{Bvh, BvhNode, BvhWorkspace};
use alloc::vec::Vec;
use smallvec::SmallVec;

const TRAVERSAL_STACK_SIZE: usize = 32;

impl Bvh {
    /*
     * Traversal of a tree against itself.
     */
    // NOTE PERF: change detection doesn’t make a huge difference in 2D (it can
    //            occasionally even make it slower!) Once we support a static/dynamic tree
    //            instead of a single tree, we might want to fully disable change detection
    //            in 2D.
    /// Traverses the Bounding Volume Test Tree of a tree against itself.
    ///
    /// The closure `f` will be called on each pair of leaf that passed the AABB intersection checks.
    /// If `CHANGE_DETECTION` is `true`, then only pairs of leaves where at least one was detected
    /// as changed during [`Self::insert_or_update_partially`] will be traversed.
    pub fn traverse_bvtt_single_tree<const CHANGE_DETECTION: bool>(
        &self,
        workspace: &mut BvhWorkspace,
        f: &mut impl FnMut(u32, u32),
    ) {
        if self.nodes.is_empty() || self.nodes[0].right.leaf_count() == 0 {
            // Not enough nodes for any overlap.
            return;
        }

        workspace.traversal_stack.clear();
        self.self_intersect_node::<CHANGE_DETECTION>(&mut workspace.traversal_stack, 0, f)
    }

    // Traverses overlaps of a single node with itself.
    // This as special case to:
    // - Ensure we don’t traverse the same branch twice.
    // - Only check the left/right overlap. Left/left and right/right checks trivially pass.
    // TODO: take change detection into account.
    fn self_intersect_node<const CHANGE_DETECTION: bool>(
        &self,
        stack: &mut Vec<u32>,
        id: u32,
        f: &mut impl FnMut(u32, u32),
    ) {
        let node = &self.nodes[id as usize];

        if CHANGE_DETECTION && !node.right.is_changed() && !node.left.is_changed() {
            return;
        }

        let left_right_intersect = node.left.intersects(&node.right);
        let left_child = node.left.children;
        let right_child = node.right.children;
        let left_is_leaf = node.left.is_leaf();
        let right_is_leaf = node.right.is_leaf();

        if (!CHANGE_DETECTION || node.left.is_changed()) && !left_is_leaf {
            self.self_intersect_node::<CHANGE_DETECTION>(stack, left_child, f);
        }

        if (!CHANGE_DETECTION || node.right.is_changed()) && !right_is_leaf {
            self.self_intersect_node::<CHANGE_DETECTION>(stack, right_child, f);
        }

        if left_right_intersect {
            match (left_is_leaf, right_is_leaf) {
                (true, true) => f(left_child, right_child),
                (true, false) => self.traverse_single_subtree::<CHANGE_DETECTION>(
                    stack,
                    &node.left,
                    right_child,
                    f,
                ),
                (false, true) => self.traverse_single_subtree::<CHANGE_DETECTION>(
                    stack,
                    &node.right,
                    left_child,
                    f,
                ),
                (false, false) => self.traverse_two_branches::<CHANGE_DETECTION>(
                    stack,
                    left_child,
                    right_child,
                    f,
                ),
            }
        }
    }

    fn traverse_two_branches<const CHANGE_DETECTION: bool>(
        &self,
        stack: &mut Vec<u32>,
        a: u32,
        b: u32,
        f: &mut impl FnMut(u32, u32),
    ) {
        let node1 = &self.nodes[a as usize];
        let node2 = &self.nodes[b as usize];

        let left1 = &node1.left;
        let right1 = &node1.right;
        let left2 = &node2.left;
        let right2 = &node2.right;

        let left_left = (!CHANGE_DETECTION || left1.is_changed() || left2.is_changed())
            && left1.intersects(left2);
        let left_right = (!CHANGE_DETECTION || left1.is_changed() || right2.is_changed())
            && left1.intersects(right2);
        let right_left = (!CHANGE_DETECTION || right1.is_changed() || left2.is_changed())
            && right1.intersects(left2);
        let right_right = (!CHANGE_DETECTION || right1.is_changed() || right2.is_changed())
            && right1.intersects(right2);

        macro_rules! dispatch(
            ($check: ident, $child_a: ident, $child_b: ident) => {
                if $check {
                    match ($child_a.is_leaf(), $child_b.is_leaf()) {
                        (true, true) => f($child_a.children, $child_b.children),
                        (true, false) => {
                            self.traverse_single_subtree::<CHANGE_DETECTION>(stack, $child_a, $child_b.children, f)
                        }
                        (false, true) => self.traverse_single_subtree::<CHANGE_DETECTION>(
                            stack,
                            $child_b,
                            $child_a.children,
                            f,
                        ),
                        (false, false) => self.traverse_two_branches::<CHANGE_DETECTION>(
                            stack,
                            $child_a.children,
                            $child_b.children,
                            f,
                        ),
                    }
                }
            }
        );

        dispatch!(left_left, left1, left2);
        dispatch!(left_right, left1, right2);
        dispatch!(right_left, right1, left2);
        dispatch!(right_right, right1, right2);
    }

    // Checks overlap between a single node and a subtree.
    fn traverse_single_subtree<const CHANGE_DETECTION: bool>(
        &self,
        stack: &mut Vec<u32>,
        node: &BvhNode,
        subtree: u32,
        f: &mut impl FnMut(u32, u32),
    ) {
        debug_assert!(stack.is_empty());

        // Since this is traversing against a single node it is more efficient to keep the leaf reference
        // around and traverse the branch using a manual stack. Left branches are traversed by the main
        // loop whereas the right branches are pushed to the stack.
        let mut curr_id = subtree;
        let node_changed = node.is_changed();

        loop {
            let curr = &self.nodes[curr_id as usize];
            let left = &curr.left;
            let right = &curr.right;
            let left_check =
                (!CHANGE_DETECTION || node_changed || left.is_changed()) && node.intersects(left);
            let right_check =
                (!CHANGE_DETECTION || node_changed || right.is_changed()) && node.intersects(right);
            let left_is_leaf = left.is_leaf();
            let right_is_leaf = right.is_leaf();
            let mut found_next = false;

            if left_check {
                if left_is_leaf {
                    f(node.children, left.children)
                } else {
                    curr_id = left.children;
                    found_next = true;
                }
            }

            if right_check {
                if right_is_leaf {
                    f(node.children, right.children)
                } else if !found_next {
                    curr_id = right.children;
                    found_next = true;
                } else {
                    // We already advanced in curr_id once, push the other
                    // branch to the stack.
                    stack.push(right.children);
                }
            }

            if !found_next {
                // Pop the stack to find the next candidate.
                if let Some(next_id) = stack.pop() {
                    curr_id = next_id;
                } else {
                    // Traversal is finished.
                    return;
                }
            }
        }
    }

    /// Parallel version of [`Self::traverse_bvtt_single_tree`], returning the reached
    /// leaf pairs instead of invoking a closure.
    ///
    /// The returned pairs are in the exact same order as the calls the sequential
    /// traversal would have made, so both versions are interchangeable, including for
    /// determinism-sensitive callers.
    #[cfg(feature = "parallel")]
    pub fn traverse_bvtt_single_tree_parallel<const CHANGE_DETECTION: bool>(
        &self,
    ) -> Vec<(u32, u32)> {
        if self.nodes.is_empty() || self.nodes[0].right.leaf_count() == 0 {
            // Not enough nodes for any overlap.
            return Vec::new();
        }

        const MAX_SPLIT_DEPTH: u32 = 6;
        self.self_intersect_node_parallel::<CHANGE_DETECTION>(0, MAX_SPLIT_DEPTH)
    }

    /// Task-parallel mirror of `self_intersect_node`: the recursions become tasks and
    /// their results are concatenated in the same order as the sequential recursion.
    #[cfg(feature = "parallel")]
    fn self_intersect_node_parallel<const CHANGE_DETECTION: bool>(
        &self,
        id: u32,
        depth: u32,
    ) -> Vec<(u32, u32)> {
        // Subtrees below this leaf count are traversed sequentially.
        const SEQ_LEAF_THRESHOLD: u32 = 512;

        let node = &self.nodes[id as usize];

        if CHANGE_DETECTION && !node.right.is_changed() && !node.left.is_changed() {
            return Vec::new();
        }

        if depth == 0 || node.left.leaf_count() + node.right.leaf_count() <= SEQ_LEAF_THRESHOLD {
            let mut out = Vec::new();
            let mut stack = Vec::new();
            self.self_intersect_node::<CHANGE_DETECTION>(&mut stack, id, &mut |a, b| {
                out.push((a, b))
            });
            return out;
        }

        let left_right_intersect = node.left.intersects(&node.right);
        let left_child = node.left.children;
        let right_child = node.right.children;
        let left_is_leaf = node.left.is_leaf();
        let right_is_leaf = node.right.is_leaf();

        let (mut left_pairs, (mut right_pairs, mut cross_pairs)) = rayon::join(
            || {
                if (!CHANGE_DETECTION || node.left.is_changed()) && !left_is_leaf {
                    self.self_intersect_node_parallel::<CHANGE_DETECTION>(left_child, depth - 1)
                } else {
                    Vec::new()
                }
            },
            || {
                rayon::join(
                    || {
                        if (!CHANGE_DETECTION || node.right.is_changed()) && !right_is_leaf {
                            self.self_intersect_node_parallel::<CHANGE_DETECTION>(
                                right_child,
                                depth - 1,
                            )
                        } else {
                            Vec::new()
                        }
                    },
                    || {
                        let mut out = Vec::new();
                        if left_right_intersect {
                            match (left_is_leaf, right_is_leaf) {
                                (true, true) => out.push((left_child, right_child)),
                                (true, false) => {
                                    let mut stack = Vec::new();
                                    self.traverse_single_subtree::<CHANGE_DETECTION>(
                                        &mut stack,
                                        &node.left,
                                        right_child,
                                        &mut |a, b| out.push((a, b)),
                                    );
                                }
                                (false, true) => {
                                    let mut stack = Vec::new();
                                    self.traverse_single_subtree::<CHANGE_DETECTION>(
                                        &mut stack,
                                        &node.right,
                                        left_child,
                                        &mut |a, b| out.push((a, b)),
                                    );
                                }
                                (false, false) => {
                                    out = self.traverse_two_branches_parallel::<CHANGE_DETECTION>(
                                        left_child,
                                        right_child,
                                        depth - 1,
                                    );
                                }
                            }
                        }
                        out
                    },
                )
            },
        );

        left_pairs.append(&mut right_pairs);
        left_pairs.append(&mut cross_pairs);
        left_pairs
    }

    /// Task-parallel mirror of `traverse_two_branches`, preserving its dispatch order
    /// (left/left, left/right, right/left, right/right).
    #[cfg(feature = "parallel")]
    fn traverse_two_branches_parallel<const CHANGE_DETECTION: bool>(
        &self,
        a: u32,
        b: u32,
        depth: u32,
    ) -> Vec<(u32, u32)> {
        const SEQ_LEAF_THRESHOLD: u32 = 512;

        let node1 = &self.nodes[a as usize];
        let node2 = &self.nodes[b as usize];
        let leaf_count = node1.left.leaf_count()
            + node1.right.leaf_count()
            + node2.left.leaf_count()
            + node2.right.leaf_count();

        if depth == 0 || leaf_count <= SEQ_LEAF_THRESHOLD {
            let mut out = Vec::new();
            let mut stack = Vec::new();
            self.traverse_two_branches::<CHANGE_DETECTION>(&mut stack, a, b, &mut |a, b| {
                out.push((a, b))
            });
            return out;
        }

        let dispatch = |child_a: &BvhNode, child_b: &BvhNode, check: bool| -> Vec<(u32, u32)> {
            let mut out = Vec::new();
            if check {
                match (child_a.is_leaf(), child_b.is_leaf()) {
                    (true, true) => out.push((child_a.children, child_b.children)),
                    (true, false) => {
                        let mut stack = Vec::new();
                        self.traverse_single_subtree::<CHANGE_DETECTION>(
                            &mut stack,
                            child_a,
                            child_b.children,
                            &mut |a, b| out.push((a, b)),
                        );
                    }
                    (false, true) => {
                        let mut stack = Vec::new();
                        self.traverse_single_subtree::<CHANGE_DETECTION>(
                            &mut stack,
                            child_b,
                            child_a.children,
                            &mut |a, b| out.push((a, b)),
                        );
                    }
                    (false, false) => {
                        out = self.traverse_two_branches_parallel::<CHANGE_DETECTION>(
                            child_a.children,
                            child_b.children,
                            depth - 1,
                        );
                    }
                }
            }
            out
        };

        let left1 = &node1.left;
        let right1 = &node1.right;
        let left2 = &node2.left;
        let right2 = &node2.right;

        let left_left = (!CHANGE_DETECTION || left1.is_changed() || left2.is_changed())
            && left1.intersects(left2);
        let left_right = (!CHANGE_DETECTION || left1.is_changed() || right2.is_changed())
            && left1.intersects(right2);
        let right_left = (!CHANGE_DETECTION || right1.is_changed() || left2.is_changed())
            && right1.intersects(left2);
        let right_right = (!CHANGE_DETECTION || right1.is_changed() || right2.is_changed())
            && right1.intersects(right2);

        let ((mut ll, mut lr), (mut rl, mut rr)) = rayon::join(
            || {
                rayon::join(
                    || dispatch(left1, left2, left_left),
                    || dispatch(left1, right2, left_right),
                )
            },
            || {
                rayon::join(
                    || dispatch(right1, left2, right_left),
                    || dispatch(right1, right2, right_right),
                )
            },
        );

        ll.append(&mut lr);
        ll.append(&mut rl);
        ll.append(&mut rr);
        ll
    }

    /// Performs a simultaneous traversal of the BVHs `self` and `other`, and yields the pairs
    /// of leaves it reached.
    ///
    /// Any node pairs failing the given `check` will be excluded from the traversal.
    pub fn leaf_pairs<'a, F: Fn(&BvhNode, &BvhNode) -> bool>(
        &'a self,
        other: &'a Self,
        check: F,
    ) -> LeafPairs<'a, F> {
        if let (Some(root1), Some(root2)) = (self.nodes.first(), other.nodes.first()) {
            let mut stack = SmallVec::default();

            if root1.left.leaf_count() > 0 && root2.right.leaf_count() > 0 {
                stack.push((&root1.left, &root2.right));
                // NOTE: we don’t need to push (&root1.left, &root2.left), it is already given as
                //       the initial value of `LeafPairs::next`.
                // stack.push((&root1.left, &root2.left))
            }

            if root1.right.leaf_count() > 0 {
                if root2.right.leaf_count() > 0 {
                    stack.push((&root1.right, &root2.right));
                }
                stack.push((&root1.right, &root2.left))
            }

            LeafPairs {
                tree1: self,
                tree2: other,
                next: Some((&root1.left, &root2.left)),
                stack,
                check,
            }
        } else {
            LeafPairs {
                tree1: self,
                tree2: other,
                next: None,
                stack: Default::default(),
                check,
            }
        }
    }
}

pub struct LeafPairs<'a, Check: Fn(&BvhNode, &BvhNode) -> bool> {
    tree1: &'a Bvh,
    tree2: &'a Bvh,
    next: Option<(&'a BvhNode, &'a BvhNode)>,
    stack: SmallVec<[(&'a BvhNode, &'a BvhNode); TRAVERSAL_STACK_SIZE]>,
    check: Check,
}

impl<'a, Check: Fn(&BvhNode, &BvhNode) -> bool> Iterator for LeafPairs<'a, Check> {
    type Item = (u32, u32);
    fn next(&mut self) -> Option<Self::Item> {
        loop {
            if self.next.is_none() {
                self.next = self.stack.pop();
            }

            let (node1, node2) = self.next.take()?;

            match (node1.is_leaf(), node2.is_leaf()) {
                (true, true) => return Some((node1.children, node2.children)),
                (true, false) => {
                    let child2 = &self.tree2.nodes[node2.children as usize];
                    if (self.check)(node1, &child2.left) {
                        self.next = Some((node1, &child2.left));
                    }
                    if (self.check)(node1, &child2.right) {
                        if self.next.is_none() {
                            self.next = Some((node1, &child2.right));
                        } else {
                            self.stack.push((node1, &child2.right));
                        }
                    }
                }
                (false, true) => {
                    let child1 = &self.tree1.nodes[node1.children as usize];
                    if (self.check)(&child1.left, node2) {
                        self.next = Some((&child1.left, node2));
                    }
                    if (self.check)(&child1.right, node2) {
                        if self.next.is_none() {
                            self.next = Some((&child1.right, node2));
                        } else {
                            self.stack.push((&child1.right, node2));
                        }
                    }
                }
                (false, false) => {
                    let child1 = &self.tree1.nodes[node1.children as usize];
                    let child2 = &self.tree2.nodes[node2.children as usize];
                    if (self.check)(&child1.left, &child2.left) {
                        self.stack.push((&child1.left, &child2.left));
                    }
                    if (self.check)(&child1.right, &child2.left) {
                        self.stack.push((&child1.right, &child2.left));
                    }
                    if (self.check)(&child1.left, &child2.right) {
                        self.stack.push((&child1.left, &child2.right));
                    }
                    if (self.check)(&child1.right, &child2.right) {
                        self.stack.push((&child1.right, &child2.right));
                    }
                }
            }
        }
    }
}
