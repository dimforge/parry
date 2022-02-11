use crate::bounding_volume::{BoundingVolume, SimdAABB, AABB};
use crate::math::Vector;
use crate::math::{Point, Real};
use crate::query::{CanonicalSplit, SplitResult};
use crate::simd::SimdReal;
use na::coordinates::X;
use simba::simd::SimdValue;

use super::utils::split_indices_wrt_dim;
use super::{IndexedData, NodeIndex, QBVHNode, QBVHProxy, QBVH};

pub struct BuilderProxies<'a, T> {
    proxies: &'a mut Vec<QBVHProxy<T>>,
    aabbs: &'a mut Vec<AABB>,
}

impl<'a, T> BuilderProxies<'a, T> {
    fn insert(&mut self, data: T, aabb: AABB)
    where
        T: IndexedData,
    {
        let index = data.index();

        if self.proxies.len() <= index {
            self.proxies.resize(index + 1, QBVHProxy::invalid());
            self.aabbs.resize(index + 1, AABB::new_invalid());
        }

        self.proxies[index] = QBVHProxy::detached(data);
        self.aabbs[index] = aabb;
    }

    fn aabbs(&self) -> &[AABB] {
        &self.aabbs
    }
}

pub trait QBVHDataSplitter<T> {
    fn split_dataset<'idx>(
        &mut self,
        subdiv_dims: [usize; 2],
        center: Point<Real>,
        indices: &'idx mut [usize],
        indices_workspace: &'idx mut Vec<usize>,
        proxies: BuilderProxies<T>,
    ) -> [&'idx mut [usize]; 4];
}

struct CenterDataSplitter;

impl<T> QBVHDataSplitter<T> for CenterDataSplitter {
    fn split_dataset<'idx>(
        &mut self,
        subdiv_dims: [usize; 2],
        center: Point<Real>,
        indices: &'idx mut [usize],
        _: &'idx mut Vec<usize>,
        proxies: BuilderProxies<T>,
    ) -> [&'idx mut [usize]; 4] {
        self.split_dataset_wo_workspace(subdiv_dims, center, indices, proxies)
    }
}

impl CenterDataSplitter {
    fn split_dataset_wo_workspace<'idx, T>(
        &mut self,
        subdiv_dims: [usize; 2],
        center: Point<Real>,
        indices: &'idx mut [usize],
        proxies: BuilderProxies<T>,
    ) -> [&'idx mut [usize]; 4] {
        // TODO: should we split wrt. the median instead of the average?
        // TODO: we should ensure each subslice contains at least 4 elements each (or less if
        // indices has less than 16 elements in the first place).
        let (left, right) = split_indices_wrt_dim(indices, &proxies.aabbs, &center, subdiv_dims[0]);

        let (left_bottom, left_top) =
            split_indices_wrt_dim(left, &proxies.aabbs, &center, subdiv_dims[1]);
        let (right_bottom, right_top) =
            split_indices_wrt_dim(right, &proxies.aabbs, &center, subdiv_dims[1]);
        [left_bottom, left_top, right_bottom, right_top]
    }
}

pub struct QbvhNonOverlappingDataSplitter<F> {
    pub canonical_split: F,
    pub epsilon: Real,
}

impl<T, F> QBVHDataSplitter<T> for QbvhNonOverlappingDataSplitter<F>
where
    T: IndexedData,
    F: FnMut(T, usize, Real, Real) -> SplitResult<(T, AABB)>,
{
    fn split_dataset<'idx>(
        &mut self,
        subdiv_dims: [usize; 2],
        center: Point<Real>,
        indices: &'idx mut [usize],
        indices_workspace: &'idx mut Vec<usize>,
        mut proxies: BuilderProxies<T>,
    ) -> [&'idx mut [usize]; 4] {
        // 1. Snap the spliting point to one fo the AABB min/max,
        // such that at least one AABB isn’t split along each dimension.
        let mut split_pt = Point::from(Vector::repeat(-Real::MAX));
        let mut split_pt_right = Point::from(Vector::repeat(Real::MAX));

        for dim in subdiv_dims {
            for i in indices.iter().copied() {
                let aabb = &proxies.aabbs[i];

                if aabb.maxs[dim] <= center[dim] && aabb.maxs[dim] > split_pt[dim] {
                    split_pt[dim] = aabb.maxs[dim];
                }

                if aabb.mins[dim] >= center[dim] && aabb.mins[dim] < split_pt_right[dim] {
                    split_pt_right[dim] = aabb.mins[dim];
                }
            }

            if (split_pt[dim] - center[dim]).abs() > (split_pt_right[dim] - center[dim]).abs() {
                split_pt[dim] = split_pt_right[dim];
            }

            if split_pt[dim] == -Real::MAX || split_pt[dim] == Real::MAX {
                // Try to at least find a splitting point that is aligned with any
                // AABB side.
                let mut candidate_min = proxies.aabbs[indices[0]].mins[dim];
                let mut candidate_max = proxies.aabbs[indices[0]].maxs[dim];
                for i in indices.iter().copied() {
                    let aabb = &proxies.aabbs[i];
                    if aabb.mins[dim] < candidate_min {
                        split_pt[dim] = candidate_min;
                        break;
                    } else if aabb.mins[dim] > candidate_min {
                        split_pt[dim] = aabb.mins[dim];
                    }

                    if aabb.maxs[dim] > candidate_max {
                        split_pt[dim] = candidate_max;
                        break;
                    } else if aabb.maxs[dim] < candidate_max {
                        split_pt[dim] = aabb.maxs[dim];
                    }
                }
            }
        }

        // If we really can’t find any splitting point along both dimensions, meaning that all the
        // aabb ranges along this dimension are equal, then split at the center.
        if (split_pt[subdiv_dims[0]] == -Real::MAX || split_pt[subdiv_dims[0]] == Real::MAX)
            && (split_pt[subdiv_dims[1]] == -Real::MAX || split_pt[subdiv_dims[1]] == Real::MAX)
        {
            split_pt = center;
        }

        // 2: Actually split the geometry.
        indices_workspace.resize(indices.len(), 0);
        indices_workspace.copy_from_slice(indices);

        for dim in subdiv_dims {
            for k in 0..indices_workspace.len() {
                let i = indices_workspace[k];
                if let SplitResult::Pair(_, _) =
                    proxies.aabbs[i].canonical_split(dim, split_pt[dim], self.epsilon)
                {
                    // The AABB was split, so we need to split the geometry too.
                    if let SplitResult::Pair((data_l, aabb_l), (data_r, aabb_r)) = (self
                        .canonical_split)(
                        proxies.proxies[i].data,
                        dim,
                        split_pt[dim],
                        self.epsilon,
                    ) {
                        indices_workspace[k] = data_l.index();
                        indices_workspace.push(data_r.index());
                        proxies.insert(data_l, aabb_l);
                        proxies.insert(data_r, aabb_r);
                    }
                }
            }
        }

        // 3: Partition the indices.
        CenterDataSplitter.split_dataset_wo_workspace(
            subdiv_dims,
            split_pt,
            indices_workspace,
            proxies,
        )
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

impl<T: IndexedData> QBVH<T> {
    /// Clears this quaternary BVH and rebuilds it from a new set of data and AABBs.
    pub fn clear_and_rebuild(
        &mut self,
        mut data_gen: impl QBVHDataGenerator<T>,
        dilation_factor: Real,
    ) {
        self.clear_and_rebuild_with_splitter(data_gen, CenterDataSplitter, dilation_factor);
    }
}

impl<T: IndexedData> QBVH<T> {
    /// Clears this quaternary BVH and rebuilds it from a new set of data and AABBs.
    pub fn clear_and_rebuild_with_splitter(
        &mut self,
        mut data_gen: impl QBVHDataGenerator<T>,
        mut splitter: impl QBVHDataSplitter<T>,
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
        let (_, aabb) = self.do_recurse_build_generic(
            &mut splitter,
            &mut indices,
            &mut aabbs,
            root_id,
            dilation_factor,
        );

        self.root_aabb = aabb;
        self.nodes[0].simd_aabb = SimdAABB::from([
            aabb,
            AABB::new_invalid(),
            AABB::new_invalid(),
            AABB::new_invalid(),
        ]);
    }

    fn do_recurse_build_generic(
        &mut self,
        splitter: &mut impl QBVHDataSplitter<T>,
        indices: &mut [usize],
        aabbs: &mut Vec<AABB>,
        parent: NodeIndex,
        dilation: Real,
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

            node.simd_aabb.dilate_by_factor(SimdReal::splat(dilation));
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

        let node = QBVHNode {
            simd_aabb: SimdAABB::new_invalid(),
            children: [0; 4], // Will be set after the recursive call
            parent,
            leaf: false,
            dirty: false,
        };

        let id = self.nodes.len() as u32;
        self.nodes.push(node);

        // Split the set along the two subdiv_dims dimensions.
        let proxies = BuilderProxies {
            proxies: &mut self.proxies,
            aabbs: aabbs,
        };

        // Recurse!
        let mut workspace = vec![]; // TODO: avoid repeated allocations?
        let splits = splitter.split_dataset(subdiv_dims, center, indices, &mut workspace, proxies);
        let n = [
            NodeIndex::new(id, 0),
            NodeIndex::new(id, 1),
            NodeIndex::new(id, 2),
            NodeIndex::new(id, 3),
        ];

        let children = [
            self.do_recurse_build_generic(splitter, splits[0], aabbs, n[0], dilation),
            self.do_recurse_build_generic(splitter, splits[1], aabbs, n[1], dilation),
            self.do_recurse_build_generic(splitter, splits[2], aabbs, n[2], dilation),
            self.do_recurse_build_generic(splitter, splits[3], aabbs, n[3], dilation),
        ];

        // Now we know the indices of the child nodes.
        self.nodes[id as usize].children =
            [children[0].0, children[1].0, children[2].0, children[3].0];
        self.nodes[id as usize].simd_aabb =
            SimdAABB::from([children[0].1, children[1].1, children[2].1, children[3].1]);
        self.nodes[id as usize]
            .simd_aabb
            .dilate_by_factor(SimdReal::splat(dilation));

        // TODO: will this chain of .merged be properly optimized?
        let my_aabb = children[0]
            .1
            .merged(&children[1].1)
            .merged(&children[2].1)
            .merged(&children[3].1);
        (id, my_aabb)
    }
}
