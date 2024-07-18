use crate::bounding_volume::{Aabb, SimdAabb};
use crate::math::Vector;
use crate::math::{Point, Real};
use crate::query::SplitResult;
use crate::simd::SimdReal;
use simba::simd::SimdValue;

use super::utils::split_indices_wrt_dim;
use super::{IndexedData, NodeIndex, Qbvh, QbvhNode, QbvhNodeFlags, QbvhProxy};

pub struct BuilderProxies<'a, LeafData> {
    proxies: &'a mut Vec<QbvhProxy<LeafData>>,
    aabbs: &'a mut Vec<Aabb>,
}

impl<'a, LeafData> BuilderProxies<'a, LeafData> {
    fn insert(&mut self, data: LeafData, aabb: Aabb)
    where
        LeafData: IndexedData,
    {
        let index = data.index();

        if self.proxies.len() <= index {
            self.proxies.resize(index + 1, QbvhProxy::invalid());
            self.aabbs.resize(index + 1, Aabb::new_invalid());
        }

        self.proxies[index] = QbvhProxy::detached(data);
        self.aabbs[index] = aabb;
    }
}

pub trait QbvhDataSplitter<LeafData> {
    fn split_dataset<'idx>(
        &mut self,
        subdiv_dims: [usize; 2],
        center: Point<Real>,
        indices: &'idx mut [usize],
        indices_workspace: &'idx mut Vec<usize>,
        proxies: BuilderProxies<LeafData>,
    ) -> [&'idx mut [usize]; 4];
}

/// A data splitter that arranges a set of Aabbs in two sets based on their center’s coordinate
/// along the split axis.
pub struct CenterDataSplitter {
    /// If all the Aabb centers have the same coordinate values along the splitting axis
    /// setting this to `true` will allow the splitter to split the Aabb set into two
    /// subsets arbitrarily.
    pub enable_fallback_split: bool,
}

impl Default for CenterDataSplitter {
    fn default() -> Self {
        Self {
            enable_fallback_split: true,
        }
    }
}

impl<LeafData> QbvhDataSplitter<LeafData> for CenterDataSplitter {
    fn split_dataset<'idx>(
        &mut self,
        subdiv_dims: [usize; 2],
        center: Point<Real>,
        indices: &'idx mut [usize],
        _: &'idx mut Vec<usize>,
        proxies: BuilderProxies<LeafData>,
    ) -> [&'idx mut [usize]; 4] {
        self.split_dataset_wo_workspace(subdiv_dims, center, indices, &*proxies.aabbs)
    }
}

impl CenterDataSplitter {
    pub(crate) fn split_dataset_wo_workspace<'idx>(
        &self,
        subdiv_dims: [usize; 2],
        center: Point<Real>,
        indices: &'idx mut [usize],
        aabbs: &[Aabb],
    ) -> [&'idx mut [usize]; 4] {
        // TODO: should we split wrt. the median instead of the average?
        // TODO: we should ensure each subslice contains at least 4 elements each (or less if
        // indices has less than 16 elements in the first place).
        let (left, right) = split_indices_wrt_dim(
            indices,
            aabbs,
            &center,
            subdiv_dims[0],
            self.enable_fallback_split,
        );

        let (left_bottom, left_top) = split_indices_wrt_dim(
            left,
            aabbs,
            &center,
            subdiv_dims[1],
            self.enable_fallback_split,
        );
        let (right_bottom, right_top) = split_indices_wrt_dim(
            right,
            aabbs,
            &center,
            subdiv_dims[1],
            self.enable_fallback_split,
        );
        [left_bottom, left_top, right_bottom, right_top]
    }
}

/// Data splitter for Qbvh construction that generates non-overlapping Aabbs at each
/// level of the tree.
///
/// This splitter assumes that no pairs of the input set of Aabb overlap (though they
/// can intersect slightly at their boundaries with an error of `epsilon`). Given this set,
/// the Qbvh constructed using this splitter will be such that no pair of intermediate nodes
/// with the same depth have overlapping Aabbs.
pub struct QbvhNonOverlappingDataSplitter<F> {
    /// The leaf data-splitting function.
    pub canonical_split: F,
    /// Allowed overlap between two leaf Aabbs.
    pub epsilon: Real,
}

impl<LeafData, F> QbvhDataSplitter<LeafData> for QbvhNonOverlappingDataSplitter<F>
where
    LeafData: IndexedData,
    F: FnMut(LeafData, usize, Real, Real, Aabb, Aabb) -> SplitResult<(LeafData, Aabb)>,
{
    fn split_dataset<'idx>(
        &mut self,
        subdiv_dims: [usize; 2],
        center: Point<Real>,
        indices: &'idx mut [usize],
        indices_workspace: &'idx mut Vec<usize>,
        mut proxies: BuilderProxies<LeafData>,
    ) -> [&'idx mut [usize]; 4] {
        // 1. Snap the splitting point to one of the Aabb min/max,
        // such that at least one Aabb isn’t split along each dimension.
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
                // Aabb side.
                let candidate_min = proxies.aabbs[indices[0]].mins[dim];
                let candidate_max = proxies.aabbs[indices[0]].maxs[dim];
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
                if let SplitResult::Pair(aabb_l, aabb_r) =
                    proxies.aabbs[i].canonical_split(dim, split_pt[dim], self.epsilon)
                {
                    // The Aabb was split, so we need to split the geometry too.
                    if let SplitResult::Pair((data_l, aabb_l), (data_r, aabb_r)) = (self
                        .canonical_split)(
                        proxies.proxies[i].data,
                        dim,
                        split_pt[dim],
                        self.epsilon,
                        aabb_l,
                        aabb_r,
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
        let center_splitter = CenterDataSplitter {
            enable_fallback_split: false,
        };

        center_splitter.split_dataset_wo_workspace(
            subdiv_dims,
            split_pt,
            indices_workspace,
            &*proxies.aabbs,
        )
    }
}

/// Trait used for generating the content of the leaves of the Qbvh acceleration structure.
pub trait QbvhDataGenerator<LeafData> {
    /// Gives an idea of the number of elements this generator contains.
    ///
    /// This is primarily used for pre-allocating some arrays for better performances.
    fn size_hint(&self) -> usize;
    /// Iterate through all the elements of this generator.
    fn for_each(&mut self, f: impl FnMut(LeafData, Aabb));
}

impl<LeafData, F> QbvhDataGenerator<LeafData> for F
where
    F: ExactSizeIterator<Item = (LeafData, Aabb)>,
{
    fn size_hint(&self) -> usize {
        self.len()
    }

    #[inline(always)]
    fn for_each(&mut self, mut f: impl FnMut(LeafData, Aabb)) {
        for (elt, aabb) in self {
            f(elt, aabb)
        }
    }
}

impl<LeafData: IndexedData> Qbvh<LeafData> {
    /// Clears this quaternary BVH and rebuilds it from a new set of data and Aabbs.
    pub fn clear_and_rebuild(
        &mut self,
        data_gen: impl QbvhDataGenerator<LeafData>,
        dilation_factor: Real,
    ) {
        self.clear_and_rebuild_with_splitter(
            data_gen,
            CenterDataSplitter::default(),
            dilation_factor,
        );
    }
}

impl<LeafData: IndexedData> Qbvh<LeafData> {
    /// Clears this quaternary BVH and rebuilds it from a new set of data and Aabbs.
    pub fn clear_and_rebuild_with_splitter(
        &mut self,
        mut data_gen: impl QbvhDataGenerator<LeafData>,
        mut splitter: impl QbvhDataSplitter<LeafData>,
        dilation_factor: Real,
    ) {
        self.free_list.clear();
        self.nodes.clear();
        self.proxies.clear();

        // Create proxies.
        let mut indices = Vec::with_capacity(data_gen.size_hint());
        let mut aabbs = vec![Aabb::new_invalid(); data_gen.size_hint()];
        self.proxies = vec![QbvhProxy::invalid(); data_gen.size_hint()];

        data_gen.for_each(|data, aabb| {
            let index = data.index();
            if index >= self.proxies.len() {
                self.proxies.resize(index + 1, QbvhProxy::invalid());
                aabbs.resize(index + 1, Aabb::new_invalid());
            }

            self.proxies[index].data = data;
            aabbs[index] = aabb;
            indices.push(index);
        });

        // Build the tree recursively.
        let root_node = QbvhNode {
            simd_aabb: SimdAabb::new_invalid(),
            children: [1, u32::MAX, u32::MAX, u32::MAX],
            parent: NodeIndex::invalid(),
            flags: QbvhNodeFlags::default(),
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
        self.nodes[0].simd_aabb = SimdAabb::from([
            aabb,
            Aabb::new_invalid(),
            Aabb::new_invalid(),
            Aabb::new_invalid(),
        ]);
    }

    fn do_recurse_build_generic(
        &mut self,
        splitter: &mut impl QbvhDataSplitter<LeafData>,
        indices: &mut [usize],
        aabbs: &mut Vec<Aabb>,
        parent: NodeIndex,
        dilation: Real,
    ) -> (u32, Aabb) {
        if indices.len() <= 4 {
            // Leaf case.
            let my_id = self.nodes.len();
            let mut leaf_aabbs = [Aabb::new_invalid(); 4];
            let mut proxy_ids = [u32::MAX; 4];

            for (k, id) in indices.iter().enumerate() {
                leaf_aabbs[k] = aabbs[*id];
                proxy_ids[k] = *id as u32;
                self.proxies[*id].node = NodeIndex::new(my_id as u32, k as u8);
            }

            let mut node = QbvhNode {
                simd_aabb: SimdAabb::from(leaf_aabbs),
                children: proxy_ids,
                parent,
                flags: QbvhNodeFlags::LEAF,
            };

            node.simd_aabb.dilate_by_factor(SimdReal::splat(dilation));
            let my_aabb = node.simd_aabb.to_merged_aabb();
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

        let center_denom = 1.0 / (indices.len() as Real);

        for i in &*indices {
            let coords = aabbs[*i].center().coords;
            center += coords * center_denom;
        }

        #[cfg(feature = "dim3")]
        {
            let variance_denom = 1.0 / ((indices.len() - 1) as Real);
            for i in &*indices {
                let dir_to_center = aabbs[*i].center() - center;
                variance += dir_to_center.component_mul(&dir_to_center) * variance_denom;
            }
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

        let node = QbvhNode {
            simd_aabb: SimdAabb::new_invalid(),
            children: [0; 4], // Will be set after the recursive call
            parent,
            flags: QbvhNodeFlags::default(),
        };

        let id = self.nodes.len() as u32;
        self.nodes.push(node);

        // Split the set along the two subdiv_dims dimensions.
        let proxies = BuilderProxies {
            proxies: &mut self.proxies,
            aabbs,
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
            SimdAabb::from([children[0].1, children[1].1, children[2].1, children[3].1]);
        self.nodes[id as usize]
            .simd_aabb
            .dilate_by_factor(SimdReal::splat(dilation));

        let my_aabb = self.nodes[id as usize].simd_aabb.to_merged_aabb();
        (id, my_aabb)
    }
}
