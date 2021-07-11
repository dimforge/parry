use crate::bounding_volume::{BoundingVolume, SimdAABB, AABB};
#[cfg(feature = "dim3")]
use crate::math::Vector;
use crate::math::{Point, Real};
use crate::partitioning::{
    SimdBestFirstVisitStatus, SimdBestFirstVisitor, SimdVisitStatus, SimdVisitor,
};
use crate::simd::{SimdReal, SIMD_WIDTH};
use crate::utils::WeightedValue;
use num::Bounded;
use simba::simd::{SimdBool, SimdValue};
use std::collections::{BinaryHeap, VecDeque};
use std::ops::Range;

/// A data to which an index is associated.
pub trait IndexedData: Copy {
    /// Creates a new default instance of `Self`.
    fn default() -> Self;
    /// Gets the index associated to `self`.
    fn index(&self) -> usize;
}

impl IndexedData for usize {
    fn default() -> Self {
        // NOTE: we use u32::MAX for compatibility
        // between 32 and 64 bit platforms.
        u32::MAX as usize
    }
    fn index(&self) -> usize {
        *self
    }
}

impl IndexedData for u32 {
    fn default() -> Self {
        u32::MAX
    }
    fn index(&self) -> usize {
        *self as usize
    }
}

impl IndexedData for u64 {
    fn default() -> Self {
        u64::MAX as u64
    }
    fn index(&self) -> usize {
        *self as usize
    }
}

/// The index of a node part of a QBVH.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
#[cfg_attr(feature = "serde-serialize", derive(Serialize, Deserialize))]
pub struct NodeIndex {
    index: u32, // Index of the addressed node in the `nodes` array.
    lane: u8,   // SIMD lane of the addressed node.
}

impl NodeIndex {
    fn new(index: u32, lane: u8) -> Self {
        Self { index, lane }
    }

    fn invalid() -> Self {
        Self {
            index: u32::MAX,
            lane: 0,
        }
    }
}

/// A SIMD node of an SIMD quad tree.
///
/// This groups four nodes of the quad-tree.
#[derive(Copy, Clone, Debug)]
#[cfg_attr(feature = "serde-serialize", derive(Serialize, Deserialize))]
struct QBVHNode {
    /// The AABBs of the qbvh nodes represented by this node.
    pub simd_aabb: SimdAABB,
    /// Index of the nodes of the 4 nodes represented by `self`.
    /// If this is a leaf, it contains the proxy ids instead.
    pub children: [u32; 4],
    /// The index of the node parent to the 4 nodes represented by `self`.
    pub parent: NodeIndex,
    /// Are the four nodes represneted by `self` leaves of the `QBVH`?
    pub leaf: bool, // TODO: pack this with the NodexIndex.lane?
    dirty: bool, // TODO: move this to a separate bitvec?
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
#[cfg_attr(feature = "serde-serialize", derive(Serialize, Deserialize))]
struct QBVHProxy<T> {
    node: NodeIndex,
    data: T, // The collider data. TODO: only set the collider generation here?
}

impl<T: IndexedData> QBVHProxy<T> {
    fn invalid() -> Self {
        Self {
            node: NodeIndex::invalid(),
            data: T::default(),
        }
    }
}

/// A quaternary bounding-volume-hierarchy.
///
/// This is a bounding-volume-hierarchy where each node has either four children or none.
#[cfg_attr(feature = "serde-serialize", derive(Serialize, Deserialize))]
#[derive(Clone, Debug)]
pub struct QBVH<T> {
    root_aabb: AABB,
    nodes: Vec<QBVHNode>,
    dirty_nodes: VecDeque<u32>,
    proxies: Vec<QBVHProxy<T>>,
}

impl<T: IndexedData> QBVH<T> {
    /// Initialize an empty quad-tree.
    pub fn new() -> Self {
        QBVH {
            root_aabb: AABB::new_invalid(),
            nodes: Vec::new(),
            dirty_nodes: VecDeque::new(),
            proxies: Vec::new(),
        }
    }

    /// Iterates mutable through all the leaf data in this QBVH.
    pub fn iter_data_mut(&mut self) -> impl Iterator<Item = (NodeIndex, &mut T)> {
        self.proxies.iter_mut().map(|p| (p.node, &mut p.data))
    }

    /// Iterate through all the leaf data in this QBVH.
    pub fn iter_data(&self) -> impl Iterator<Item = (NodeIndex, &T)> {
        self.proxies.iter().map(|p| (p.node, &p.data))
    }

    /// The AABB of the root of this tree.
    pub fn root_aabb(&self) -> &AABB {
        &self.root_aabb
    }

    /// Returns the data associated to a given leaf.
    ///
    /// Returns `None` if the provided node ID does not identify a leaf.
    pub fn leaf_data(&mut self, node_id: NodeIndex) -> Option<T> {
        let node = self.nodes.get(node_id.index as usize)?;

        if !node.leaf {
            return None;
        }

        let proxy = self
            .proxies
            .get(node.children[node_id.lane as usize] as usize)?;
        Some(proxy.data)
    }

    /// Clears this quad-tree and rebuilds it from a new set of data and AABBs.
    pub fn clear_and_rebuild(
        &mut self,
        mut data_gen: impl QBVHDataGenerator<T>,
        dilation_factor: Real,
    ) {
        self.nodes.clear();
        self.proxies.clear();

        // Create proxies.
        let mut indices = Vec::with_capacity(data_gen.size_hint());
        let mut aabbs = vec![AABB::new_invalid(); data_gen.size_hint()];
        self.proxies = vec![QBVHProxy::invalid(); data_gen.size_hint()];

        data_gen.for_each(|data, aabb| {
            let index = data.index();
            if index >= self.proxies.len() {
                self.proxies.resize(index + 1, QBVHProxy::invalid());
                aabbs.resize(index + 1, AABB::new_invalid());
            }

            self.proxies[index].data = data;
            aabbs[index] = aabb;
            indices.push(index);
        });

        // Build the tree recursively.
        let root_node = QBVHNode {
            simd_aabb: SimdAABB::new_invalid(),
            children: [1, u32::MAX, u32::MAX, u32::MAX],
            parent: NodeIndex::invalid(),
            leaf: false,
            dirty: false,
        };

        self.nodes.push(root_node);
        let root_id = NodeIndex::new(0, 0);
        let (_, aabb) = self.do_recurse_build(&mut indices, &aabbs, root_id, dilation_factor);
        self.root_aabb = aabb;
        self.nodes[0].simd_aabb = SimdAABB::from([
            aabb,
            AABB::new_invalid(),
            AABB::new_invalid(),
            AABB::new_invalid(),
        ]);
    }

    /// Marks a piece of data as dirty so it can be updated during the next
    /// call to `self.update`.
    pub fn pre_update(&mut self, data: T) {
        let id = data.index();
        let node_id = self.proxies[id].node.index;
        let node = &mut self.nodes[node_id as usize];
        if !node.dirty {
            node.dirty = true;
            self.dirty_nodes.push_back(node_id);
        }
    }

    /// Update all the nodes that have been marked as dirty by `self.pre_update`.
    pub fn update<F>(&mut self, aabb_builder: F, dilation_factor: Real)
    where
        F: Fn(&T) -> AABB,
    {
        // Loop on the dirty leaves.
        let dilation_factor = SimdReal::splat(dilation_factor);

        while let Some(id) = self.dirty_nodes.pop_front() {
            // NOTE: this will data the case where we reach the root of the tree.
            if let Some(node) = self.nodes.get(id as usize) {
                // Compute the new aabb.
                let mut new_aabbs = [AABB::new_invalid(); SIMD_WIDTH];
                for (child_id, new_aabb) in node.children.iter().zip(new_aabbs.iter_mut()) {
                    if node.leaf {
                        // We are in a leaf: compute the AABBs.
                        if let Some(proxy) = self.proxies.get(*child_id as usize) {
                            *new_aabb = aabb_builder(&proxy.data);
                        }
                    } else {
                        // We are in an internal node: compute the children's AABBs.
                        if let Some(node) = self.nodes.get(*child_id as usize) {
                            *new_aabb = node.simd_aabb.to_merged_aabb();
                        }
                    }
                }

                let node = &mut self.nodes[id as usize];
                let new_simd_aabb = SimdAABB::from(new_aabbs);
                if !node.simd_aabb.contains(&new_simd_aabb).all() {
                    node.simd_aabb = new_simd_aabb;
                    node.simd_aabb.dilate_by_factor(dilation_factor);
                    self.dirty_nodes.push_back(node.parent.index);
                }
                node.dirty = false;
            }
        }
    }

    fn do_recurse_build(
        &mut self,
        indices: &mut [usize],
        aabbs: &[AABB],
        parent: NodeIndex,
        dilation_factor: Real,
    ) -> (u32, AABB) {
        if indices.len() <= 4 {
            // Leaf case.
            let my_id = self.nodes.len();
            let mut my_aabb = AABB::new_invalid();
            let mut leaf_aabbs = [AABB::new_invalid(); 4];
            let mut proxy_ids = [u32::MAX; 4];

            for (k, id) in indices.iter().enumerate() {
                my_aabb.merge(&aabbs[*id]);
                leaf_aabbs[k] = aabbs[*id];
                proxy_ids[k] = *id as u32;
                self.proxies[*id].node = NodeIndex::new(my_id as u32, k as u8);
            }

            let mut node = QBVHNode {
                simd_aabb: SimdAABB::from(leaf_aabbs),
                children: proxy_ids,
                parent,
                leaf: true,
                dirty: false,
            };

            node.simd_aabb
                .dilate_by_factor(SimdReal::splat(dilation_factor));
            self.nodes.push(node);
            return (my_id as u32, my_aabb);
        }

        // Compute the center and variance along each dimension.
        // In 3D we compute the variance to not-subdivide the dimension with lowest variance.
        // Therefore variance computation is not needed in 2D because we only have 2 dimension
        // to split in the first place.
        let mut center = Point::origin();
        #[cfg(feature = "dim3")]
        let mut variance = Vector::zeros();

        let denom = 1.0 / (indices.len() as Real);

        for i in &*indices {
            let coords = aabbs[*i].center().coords;
            center += coords * denom;
            #[cfg(feature = "dim3")]
            {
                variance += coords.component_mul(&coords) * denom;
            }
        }

        #[cfg(feature = "dim3")]
        {
            variance = variance - center.coords.component_mul(&center.coords);
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

        // Split the set along the two subdiv_dims dimensions.
        // TODO: should we split wrt. the median instead of the average?
        // TODO: we should ensure each subslice contains at least 4 elements each (or less if
        // indices has less than 16 elements in the first place.
        let (left, right) = split_indices_wrt_dim(indices, &aabbs, &center, subdiv_dims[0]);

        let (left_bottom, left_top) = split_indices_wrt_dim(left, &aabbs, &center, subdiv_dims[1]);
        let (right_bottom, right_top) =
            split_indices_wrt_dim(right, &aabbs, &center, subdiv_dims[1]);

        // println!(
        //     "Recursing on children: {}, {}, {}, {}",
        //     left_bottom.len(),
        //     left_top.len(),
        //     right_bottom.len(),
        //     right_top.len()
        // );

        let node = QBVHNode {
            simd_aabb: SimdAABB::new_invalid(),
            children: [0; 4], // Will be set after the recursive call
            parent,
            leaf: false,
            dirty: false,
        };

        let id = self.nodes.len() as u32;
        self.nodes.push(node);

        // Recurse!
        let a = self.do_recurse_build(left_bottom, aabbs, NodeIndex::new(id, 0), dilation_factor);
        let b = self.do_recurse_build(left_top, aabbs, NodeIndex::new(id, 1), dilation_factor);
        let c = self.do_recurse_build(right_bottom, aabbs, NodeIndex::new(id, 2), dilation_factor);
        let d = self.do_recurse_build(right_top, aabbs, NodeIndex::new(id, 3), dilation_factor);

        // Now we know the indices of the grand-nodes.
        self.nodes[id as usize].children = [a.0, b.0, c.0, d.0];
        self.nodes[id as usize].simd_aabb = SimdAABB::from([a.1, b.1, c.1, d.1]);
        self.nodes[id as usize]
            .simd_aabb
            .dilate_by_factor(SimdReal::splat(dilation_factor));

        // TODO: will this chain of .merged be properly optimized?
        let my_aabb = a.1.merged(&b.1).merged(&c.1).merged(&d.1);
        (id, my_aabb)
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
}

/// Trait used for generating the content of the leaves of the QBVH acceleration structure.
pub trait QBVHDataGenerator<T> {
    /// Gives an idea of the number of elements this generator contains.
    ///
    /// This is primarily used for pre-allocating some arrays for better performances.
    fn size_hint(&self) -> usize;
    /// Iterate through all the elements of this generator.
    fn for_each(&mut self, f: impl FnMut(T, AABB));
}

impl<T, F> QBVHDataGenerator<T> for F
where
    F: ExactSizeIterator<Item = (T, AABB)>,
{
    fn size_hint(&self) -> usize {
        self.len()
    }

    #[inline(always)]
    fn for_each(&mut self, mut f: impl FnMut(T, AABB)) {
        for (elt, aabb) in self {
            f(elt, aabb)
        }
    }
}

#[allow(dead_code)]
struct QBVHIncrementalBuilderStep {
    range: Range<usize>,
    parent: NodeIndex,
}

#[allow(dead_code)]
struct QBVHIncrementalBuilder<T> {
    qbvh: QBVH<T>,
    to_insert: Vec<QBVHIncrementalBuilderStep>,
    aabbs: Vec<AABB>,
    indices: Vec<usize>,
}

#[allow(dead_code)]
impl<T: IndexedData> QBVHIncrementalBuilder<T> {
    pub fn new() -> Self {
        Self {
            qbvh: QBVH::new(),
            to_insert: Vec::new(),
            aabbs: Vec::new(),
            indices: Vec::new(),
        }
    }

    pub fn update_single_depth(&mut self) {
        if let Some(to_insert) = self.to_insert.pop() {
            let indices = &mut self.indices[to_insert.range];

            // Leaf case.
            if indices.len() <= 4 {
                let id = self.qbvh.nodes.len();
                let mut aabb = AABB::new_invalid();
                let mut leaf_aabbs = [AABB::new_invalid(); 4];
                let mut proxy_ids = [u32::MAX; 4];

                for (k, id) in indices.iter().enumerate() {
                    aabb.merge(&self.aabbs[*id]);
                    leaf_aabbs[k] = self.aabbs[*id];
                    proxy_ids[k] = *id as u32;
                }

                let node = QBVHNode {
                    simd_aabb: SimdAABB::from(leaf_aabbs),
                    children: proxy_ids,
                    parent: to_insert.parent,
                    leaf: true,
                    dirty: false,
                };

                self.qbvh.nodes[to_insert.parent.index as usize].children
                    [to_insert.parent.lane as usize] = id as u32;
                self.qbvh.nodes[to_insert.parent.index as usize]
                    .simd_aabb
                    .replace(to_insert.parent.lane as usize, aabb);
                self.qbvh.nodes.push(node);
                return;
            }

            // Compute the center and variance along each dimension.
            // In 3D we compute the variance to not-subdivide the dimension with lowest variance.
            // Therefore variance computation is not needed in 2D because we only have 2 dimension
            // to split in the first place.
            let mut center = Point::origin();
            #[cfg(feature = "dim3")]
            let mut variance = Vector::zeros();

            let denom = 1.0 / (indices.len() as Real);
            let mut aabb = AABB::new_invalid();

            for i in &*indices {
                let coords = self.aabbs[*i].center().coords;
                aabb.merge(&self.aabbs[*i]);
                center += coords * denom;
                #[cfg(feature = "dim3")]
                {
                    variance += coords.component_mul(&coords) * denom;
                }
            }

            #[cfg(feature = "dim3")]
            {
                variance = variance - center.coords.component_mul(&center.coords);
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

            // Split the set along the two subdiv_dims dimensions.
            // TODO: should we split wrt. the median instead of the average?
            // TODO: we should ensure each subslice contains at least 4 elements each (or less if
            // indices has less than 16 elements in the first place.
            let (left, right) =
                split_indices_wrt_dim(indices, &self.aabbs, &center, subdiv_dims[0]);

            let (left_bottom, left_top) =
                split_indices_wrt_dim(left, &self.aabbs, &center, subdiv_dims[1]);
            let (right_bottom, right_top) =
                split_indices_wrt_dim(right, &self.aabbs, &center, subdiv_dims[1]);

            let node = QBVHNode {
                simd_aabb: SimdAABB::new_invalid(),
                children: [0; 4], // Will be set after the recursive call
                parent: to_insert.parent,
                leaf: false,
                dirty: false,
            };

            let id = self.qbvh.nodes.len() as u32;
            self.qbvh.nodes.push(node);

            // Recurse!
            let a = left_bottom.len();
            let b = a + left_top.len();
            let c = b + right_bottom.len();
            let d = c + right_top.len();
            self.to_insert.push(QBVHIncrementalBuilderStep {
                range: 0..a,
                parent: NodeIndex::new(id, 0),
            });
            self.to_insert.push(QBVHIncrementalBuilderStep {
                range: a..b,
                parent: NodeIndex::new(id, 1),
            });
            self.to_insert.push(QBVHIncrementalBuilderStep {
                range: b..c,
                parent: NodeIndex::new(id, 2),
            });
            self.to_insert.push(QBVHIncrementalBuilderStep {
                range: c..d,
                parent: NodeIndex::new(id, 3),
            });

            self.qbvh.nodes[to_insert.parent.index as usize].children
                [to_insert.parent.lane as usize] = id as u32;
            self.qbvh.nodes[to_insert.parent.index as usize]
                .simd_aabb
                .replace(to_insert.parent.lane as usize, aabb);
        }
    }
}

fn split_indices_wrt_dim<'a>(
    indices: &'a mut [usize],
    aabbs: &[AABB],
    split_point: &Point<Real>,
    dim: usize,
) -> (&'a mut [usize], &'a mut [usize]) {
    let mut icurr = 0;
    let mut ilast = indices.len();

    // The loop condition we can just do 0..indices.len()
    // instead of the test icurr < ilast because we know
    // we will iterate exactly once per index.
    for _ in 0..indices.len() {
        let i = indices[icurr];
        let center = aabbs[i].center();

        if center[dim] > split_point[dim] {
            ilast -= 1;
            indices.swap(icurr, ilast);
        } else {
            icurr += 1;
        }
    }

    if icurr == 0 || icurr == indices.len() {
        // We don't want to return one empty set. But
        // this can happen if all the coordinates along the
        // given dimension are equal.
        // In this is the case, we just split in the middle.
        let half = indices.len() / 2;
        indices.split_at_mut(half)
    } else {
        indices.split_at_mut(icurr)
    }
}

#[cfg(test)]
mod test {
    use crate::bounding_volume::AABB;
    use crate::math::{Point, Vector};
    use crate::partitioning::QBVH;

    #[test]
    fn multiple_identical_aabb_stack_overflow() {
        // A stack overflow was caused during the construction of the
        // QBVH with more than four AABB with the same center.
        let aabb = AABB::new(Point::origin(), Vector::repeat(1.0).into());

        for k in 0u32..20 {
            let mut tree = QBVH::new();
            tree.clear_and_rebuild((0..k).map(|i| (i, aabb)), 0.0);
        }
    }
}
