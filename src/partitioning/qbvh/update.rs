use crate::bounding_volume::{BoundingVolume, SimdAABB, AABB};
#[cfg(feature = "dim3")]
use crate::math::Vector;
use crate::math::{Point, Real};
use crate::simd::{SimdReal, SIMD_WIDTH};
use simba::simd::{SimdBool, SimdValue};
use std::ops::Range;

use super::{utils::split_indices_wrt_dim, IndexedData, NodeIndex, QBVHNode, QBVH};

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
                split_indices_wrt_dim(indices, &self.aabbs, &center, subdiv_dims[0], true);

            let (left_bottom, left_top) =
                split_indices_wrt_dim(left, &self.aabbs, &center, subdiv_dims[1], true);
            let (right_bottom, right_top) =
                split_indices_wrt_dim(right, &self.aabbs, &center, subdiv_dims[1], true);

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

impl<T: IndexedData> QBVH<T> {
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
            // NOTE: this will deal with the case where we reach the root of the tree.
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
}
