use crate::bounding_volume::{BoundingVolume, SimdAABB, AABB};
#[cfg(feature = "dim3")]
use crate::math::Vector;
use crate::math::{Point, Real};
use crate::simd::SimdReal;
use either::Either;
use simba::simd::SimdValue;

use super::utils::split_indices_wrt_dim;
use super::{IndexedData, NodeIndex, QBVHNode, QBVHProxy, QBVH};

pub struct BuilderProxies<'a, T> {
    proxies: &'a mut Vec<QBVHProxy<T>>,
    aabbs: &'a mut Vec<AABB>,
}

impl<'a, T> BuilderProxies<'a, T> {
    #[must_use]
    #[inline(always)]
    fn get(&self, i: usize) -> (&T, &AABB) {
        (&self.proxies[i].data, &self.aabbs[i])
    }

    fn insert(&mut self, data: T, aabb: AABB) {
        self.proxies.push(QBVHProxy::detached(data));
        self.aabbs.push(aabb);
    }
}

pub trait QBVHDataSplitter<T> {
    fn split_dataset(
        &mut self,
        indices: &mut [usize],
        proxies: &mut BuilderProxies<T>,
    ) -> Either<[&mut [usize]; 4], [Vec<usize>; 4]>;
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
}
