use crate::bounding_volume::Aabb;
use crate::math::{Point, Real, Vector, DIM};
use crate::shape::voxels::voxels_chunk::{VoxelsChunk, VoxelsChunkHeader};
use crate::shape::{VoxelState, Voxels};
use crate::utils::hashmap::Entry;
use alloc::vec;

impl Voxels {
    /// Sets the size of each voxel along each local coordinate axis.
    ///
    /// Since the internal spatial acceleration structure needs to be updated, this
    /// operation runs in `O(n)` time, where `n` is the number of voxels.
    pub fn set_voxel_size(&mut self, new_size: Vector<Real>) {
        let scale = new_size.component_div(&self.voxel_size);
        self.chunk_bvh.scale(&scale);
        self.voxel_size = new_size;
    }

    /// Inserts or remove a voxel from this shape.
    ///
    /// Return the previous `VoxelState` of this voxel.
    pub fn set_voxel(&mut self, key: Point<i32>, is_filled: bool) -> VoxelState {
        let (chunk_key, id_in_chunk) = Self::chunk_key_and_id_in_chunk(key);
        let header_entry = self.chunk_headers.entry(chunk_key);

        if !is_filled && matches!(header_entry, Entry::Vacant(_)) {
            // The voxel is already empty (it doesn’t exist at all).
            // Nothing more to do.
            return VoxelState::EMPTY;
        }

        let chunk_header = header_entry.or_insert_with(|| {
            let id = self.free_chunks.pop().unwrap_or_else(|| {
                self.chunks.push(VoxelsChunk::default());
                self.chunk_keys.push(chunk_key);
                self.chunks.len() - 1
            });

            self.chunk_keys[id] = chunk_key;
            self.chunk_bvh
                .insert(VoxelsChunk::aabb(&chunk_key, &self.voxel_size), id as u32);
            VoxelsChunkHeader { id, len: 0 }
        });
        let chunk_id = chunk_header.id;

        let prev = self.chunks[chunk_id].states[id_in_chunk];
        let new_is_empty = !is_filled;

        if prev.is_empty() ^ new_is_empty {
            let can_remove_chunk = if new_is_empty {
                chunk_header.len -= 1;
                chunk_header.len == 0
            } else {
                chunk_header.len += 1;
                false
            };

            self.chunks[chunk_id].states[id_in_chunk] =
                self.update_neighbors_state(key, new_is_empty);

            if can_remove_chunk {
                self.chunk_bvh.remove(chunk_id as u32);

                #[cfg(feature = "enhanced-determinism")]
                let _ = self.chunk_headers.swap_remove(&chunk_key);
                #[cfg(not(feature = "enhanced-determinism"))]
                let _ = self.chunk_headers.remove(&chunk_key);

                self.free_chunks.push(chunk_id);
                self.chunk_keys[chunk_id] = VoxelsChunk::INVALID_CHUNK_KEY;
            }
        }

        prev
    }

    /// Crops in-place the voxel shape with a rectangular domain.
    ///
    /// This removes every voxels out of the `[domain_mins, domain_maxs]` bounds.
    pub fn crop(&mut self, domain_mins: Point<i32>, domain_maxs: Point<i32>) {
        // TODO PERF: this could be done more efficiently.
        if let Some(new_shape) = self.cropped(domain_mins, domain_maxs) {
            *self = new_shape;
        }
    }

    /// Returns a cropped version of this voxel shape with a rectangular domain.
    ///
    /// This removes every voxels out of the `[domain_mins, domain_maxs]` bounds.
    pub fn cropped(&self, domain_mins: Point<i32>, domain_maxs: Point<i32>) -> Option<Self> {
        // TODO PERF: can be optimized significantly.
        let mut in_box = vec![];
        for vox in self.voxels() {
            if !vox.state.is_empty()
                && grid_aabb_contains_point(&domain_mins, &domain_maxs, &vox.grid_coords)
            {
                in_box.push(vox.grid_coords);
            }
        }

        if !in_box.is_empty() {
            Some(Voxels::new(self.voxel_size, &in_box))
        } else {
            None
        }
    }

    /// Splits this voxels shape into two subshapes.
    ///
    /// The first subshape contains all the voxels which centers are inside the `aabb`.
    /// The second subshape contains all the remaining voxels.
    pub fn split_with_box(&self, aabb: &Aabb) -> (Option<Self>, Option<Self>) {
        // TODO PERF: can be optimized significantly.
        let mut in_box = vec![];
        let mut rest = vec![];
        for vox in self.voxels() {
            if !vox.state.is_empty() {
                if aabb.contains_local_point(&vox.center) {
                    in_box.push(vox.grid_coords);
                } else {
                    rest.push(vox.grid_coords);
                }
            }
        }

        let in_box = if !in_box.is_empty() {
            Some(Voxels::new(self.voxel_size, &in_box))
        } else {
            None
        };

        let rest = if !rest.is_empty() {
            Some(Voxels::new(self.voxel_size, &rest))
        } else {
            None
        };

        (in_box, rest)
    }
}

fn grid_aabb_contains_point(mins: &Point<i32>, maxs: &Point<i32>, point: &Point<i32>) -> bool {
    for i in 0..DIM {
        if point[i] < mins[i] || point[i] > maxs[i] {
            return false;
        }
    }

    true
}
