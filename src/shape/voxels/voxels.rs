use super::{
    VoxelIndex, EMPTY_FACE_MASK, FACES_TO_FEATURE_MASKS, FACES_TO_OCTANT_MASKS,
    FACES_TO_VOXEL_TYPES, INTERIOR_FACE_MASK,
};
use crate::bounding_volume::{Aabb, BoundingVolume};
use crate::math::{Point, Real, Vector};
use crate::partitioning::{Bvh, BvhBuildStrategy, BvhNode};
use crate::shape::voxels::voxels_chunk::{VoxelsChunk, VoxelsChunkHeader};
use crate::shape::VoxelsChunkRef;
use crate::utils::hashmap::HashMap;
use alloc::{vec, vec::Vec};
#[cfg(not(feature = "std"))]
use na::ComplexField;

/// Categorization of a voxel based on its neighbors.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum VoxelType {
    /// The voxel is empty.
    Empty,
    /// The voxel is a vertex if all three coordinate axis directions have at
    /// least one empty neighbor.
    Vertex,
    /// The voxel is on an edge if it has non-empty neighbors in both directions of
    /// a single coordinate axis.
    #[cfg(feature = "dim3")]
    Edge,
    /// The voxel is on an edge if it has non-empty neighbors in both directions of
    /// two coordinate axes.
    Face,
    /// The voxel is on an edge if it has non-empty neighbors in both directions of
    /// all coordinate axes.
    Interior,
}

#[derive(Clone, Copy, Debug, Default, Eq, Hash, Ord, PartialEq, PartialOrd)]
/// The status of the cell of an heightfield.
pub struct AxisMask(u8);

bitflags::bitflags! {
    /// Flags for identifying signed directions along coordinate axes, or faces of a voxel.
    impl AxisMask: u8 {
        /// The direction or face along the `+x` coordinate axis.
        const X_POS = 1 << 0;
        /// The direction or face along the `-x` coordinate axis.
        const X_NEG = 1 << 1;
        /// The direction or face along the `+y` coordinate axis.
        const Y_POS = 1 << 2;
        /// The direction or face along the `-y` coordinate axis.
        const Y_NEG = 1 << 3;
        /// The direction or face along the `+z` coordinate axis.
        #[cfg(feature= "dim3")]
        const Z_POS = 1 << 4;
        /// The direction or face along the `-z` coordinate axis.
        #[cfg(feature= "dim3")]
        const Z_NEG = 1 << 5;
    }
}

/// Indicates the local shape of a voxel on each octant.
///
/// This provides geometric information of the shape’s exposed features on each octant.
// This is an alternative to `FACES_TO_FEATURE_MASKS` that can be more convenient for some
// collision-detection algorithms.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct OctantPattern;

// NOTE: it is important that the max value of any OctantPattern variant
//       is 7 because we don’t allocate more than 3 bits to store it in
//      `FACES_TO_OCTANT_MASKS`.
/// Indicates the local shape of a voxel on each octant.
///
/// This provides geometric information of the shape’s exposed features on each octant.
// This is an alternative to `FACES_TO_FEATURE_MASKS` that can be more convenient for some
// collision-detection algorithms.
#[cfg(feature = "dim3")]
impl OctantPattern {
    /// The voxel doesn't have any exposed feature on the octant with this mask.
    pub const INTERIOR: u32 = 0;
    /// The voxel has an exposed vertex on the octant with this mask.
    pub const VERTEX: u32 = 1;
    /// The voxel has an exposed edges with direction X on the octant with this mask.
    pub const EDGE_X: u32 = 2;
    /// The voxel has an exposed edges with direction Y on the octant with this mask.
    pub const EDGE_Y: u32 = 3;
    /// The voxel has an exposed edges with direction Z on the octant with this mask.
    pub const EDGE_Z: u32 = 4;
    /// The voxel has an exposed face with normal X on the octant with this mask.
    pub const FACE_X: u32 = 5;
    /// The voxel has an exposed face with normal Y on the octant with this mask.
    pub const FACE_Y: u32 = 6;
    /// The voxel has an exposed face with normal Z on the octant with this mask.
    pub const FACE_Z: u32 = 7;
}

// NOTE: it is important that the max value of any OctantPattern variant
//       is 7 because we don’t allocate more than 3 bits to store it in
//      `FACES_TO_OCTANT_MASKS`.
/// Indicates the local shape of a voxel on each octant.
///
/// This provides geometric information of the shape’s exposed features on each octant.
/// This is an alternative to `FACES_TO_FEATURE_MASKS` that can be more convenient for some
/// collision-detection algorithms.
#[cfg(feature = "dim2")]
impl OctantPattern {
    /// The voxel doesn't have any exposed feature on the octant with this mask.
    pub const INTERIOR: u32 = 0;
    /// The voxel has an exposed vertex on the octant with this mask.
    pub const VERTEX: u32 = 1;
    /// The voxel has an exposed face with normal X on the octant with this mask.
    pub const FACE_X: u32 = 2;
    /// The voxel has an exposed face with normal Y on the octant with this mask.
    pub const FACE_Y: u32 = 3;
}

// The local neighborhood information is encoded in a 8-bits number in groups of two bits
// per coordinate axis: `0bwwzzyyxx`. In each group of two bits, e.g. `xx`, the rightmost (resp.
// leftmost) bit set to 1 means that the neighbor voxel with coordinate `+1` (resp `-1`) relative
// to the current voxel along the `x` axis is filled. If the bit is 0, then the corresponding
// neighbor is empty. See the `AxisMask` bitflags.
// For example, in 2D, the mask `0b00_00_10_01` matches the following configuration (assuming +y goes
// up, and +x goes right):
//
// ```txt
//  0 0 0
//  0 x 1
//  0 1 0
// ```
//
// The special value `0b01000000` indicates that the voxel is empty.
// And the value `0b00111111` (`0b00001111` in 2D) indicates that the voxel is an interior voxel (its whole neighborhood
// is filled).
/// A description of the local neighborhood of a voxel.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
#[cfg_attr(feature = "serde-serialize", derive(Serialize, Deserialize))]
pub struct VoxelState(pub(super) u8);

impl VoxelState {
    /// The value of empty voxels.
    pub const EMPTY: VoxelState = VoxelState(EMPTY_FACE_MASK);
    /// The value of a voxel with non-empty neighbors in all directions.
    pub const INTERIOR: VoxelState = VoxelState(INTERIOR_FACE_MASK);

    pub(crate) const fn new(state: u8) -> Self {
        Self(state)
    }

    /// Is this voxel empty?
    pub const fn is_empty(self) -> bool {
        self.0 == EMPTY_FACE_MASK
    }

    /// A bit mask indicating which faces of the voxel don’t have any
    /// adjacent non-empty voxel.
    pub const fn free_faces(self) -> AxisMask {
        if self.0 == INTERIOR_FACE_MASK || self.0 == EMPTY_FACE_MASK {
            AxisMask::empty()
        } else {
            AxisMask::from_bits_truncate((!self.0) & INTERIOR_FACE_MASK)
        }
    }

    /// The [`VoxelType`] of this voxel.
    pub const fn voxel_type(self) -> VoxelType {
        FACES_TO_VOXEL_TYPES[self.0 as usize]
    }

    // Bitmask indicating what vertices, edges, or faces of the voxel are "free".
    pub(crate) const fn feature_mask(self) -> u16 {
        FACES_TO_FEATURE_MASKS[self.0 as usize]
    }

    pub(crate) const fn octant_mask(self) -> u32 {
        FACES_TO_OCTANT_MASKS[self.0 as usize]
    }
}

/// Information associated to a voxel.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct VoxelData {
    /// The temporary index in the internal voxels’ storage.
    ///
    /// This index can be invalidated after a call to [`Voxels::set_voxel`], or
    /// [`Voxels::crop`].
    pub linear_id: VoxelIndex,
    /// The voxel’s integer grid coordinates.
    pub grid_coords: Point<i32>,
    /// The voxel’s center position in the local-space of the [`Voxels`] shape it is part of.
    pub center: Point<Real>,
    /// The voxel’s state, indicating if it’s empty or full.
    pub state: VoxelState,
}

/// A shape made of axis-aligned, uniformly sized, cubes (aka. voxels).
///
/// This shape is specialized to handle voxel worlds and voxelized obojects efficiently why ensuring
/// that collision-detection isn’t affected by the so-called "internal edges problem" that can create
/// artifacts when another object rolls or slides against a flat voxelized surface.
///
/// The internal storage is compact (but not sparse at the moment), storing only one byte per voxel
/// in the allowed domain. This has a generally smaller memory footprint than a mesh representation
/// of the voxels.
#[derive(Clone, Debug)]
#[cfg_attr(feature = "serde-serialize", derive(Serialize, Deserialize))]
pub struct Voxels {
    /// A BVH of chunk keys.
    ///
    /// The bounding boxes are the ones of the chunk’s voxels **keys**. This is equivalent to a bvh
    /// of the chunks with a uniform voxel size of 1.
    pub(super) chunk_bvh: Bvh,
    pub(super) chunk_headers: HashMap<Point<i32>, VoxelsChunkHeader>,
    pub(super) chunk_keys: Vec<Point<i32>>,
    pub(super) chunks: Vec<VoxelsChunk>,
    pub(super) free_chunks: Vec<usize>,
    pub(super) voxel_size: Vector<Real>,
}

impl Voxels {
    /// Initializes a voxel shapes from the voxels grid coordinates.
    ///
    /// Each voxel will have its bottom-left-back corner located at
    /// `grid_coordinates * voxel_size`; and its center at `(grid_coordinates + 0.5) * voxel_size`.
    pub fn new(voxel_size: Vector<Real>, grid_coordinates: &[Point<i32>]) -> Self {
        let mut result = Self {
            chunk_bvh: Bvh::new(),
            chunk_headers: HashMap::default(),
            chunk_keys: vec![],
            chunks: vec![],
            free_chunks: vec![],
            voxel_size,
        };

        for vox in grid_coordinates {
            let (chunk_key, id_in_chunk) = Self::chunk_key_and_id_in_chunk(*vox);
            let chunk_header = result.chunk_headers.entry(chunk_key).or_insert_with(|| {
                let id = result.chunks.len();
                result.chunks.push(VoxelsChunk::default());
                result.chunk_keys.push(chunk_key);
                VoxelsChunkHeader { id, len: 0 }
            });
            chunk_header.len += 1;
            result.chunks[chunk_header.id].states[id_in_chunk] = VoxelState::INTERIOR;
        }

        result.chunk_bvh = Bvh::from_iter(
            BvhBuildStrategy::Ploc,
            result.chunk_headers.iter().map(|(chunk_key, chunk_id)| {
                (
                    chunk_id.id,
                    VoxelsChunk::aabb(chunk_key, &result.voxel_size),
                )
            }),
        );

        result.recompute_all_voxels_states();
        result
    }

    /// Computes a voxels shape from the set of `points`.
    ///
    /// The points are mapped to a regular grid centered at the provided point with smallest
    /// coordinates, and with grid cell size equal to `scale`. It is OK if multiple points
    /// fall into the same grid cell.
    pub fn from_points(voxel_size: Vector<Real>, points: &[Point<Real>]) -> Self {
        let voxels: Vec<_> = points
            .iter()
            .map(|pt| {
                Point::from(
                    pt.coords
                        .component_div(&voxel_size)
                        .map(|x| x.floor() as i32),
                )
            })
            .collect();
        Self::new(voxel_size, &voxels)
    }

    pub(crate) fn chunk_bvh(&self) -> &Bvh {
        &self.chunk_bvh
    }

    /// The semi-open range of voxels in shape.
    ///
    /// This provides conservative bounds on the range of voxel indices that might be set to filled.
    pub fn domain(&self) -> [Point<i32>; 2] {
        let aabb = self.chunk_bvh.root_aabb();

        // NOTE that we shift the AABB’s bounds so the endpoint matches a voxel center
        //      to avoid rounding errors.
        let half_sz = self.voxel_size() / 2.0;
        let mins = self.voxel_at_point(aabb.mins + half_sz);
        // + 1 because the range is semi-open.
        let maxs = self.voxel_at_point(aabb.maxs - half_sz) + Vector::repeat(1);
        [mins, maxs]
    }

    // /// The number of voxels along each coordinate axis.
    // pub fn dimensions(&self) -> Vector<u32> {
    //     (self.domain_maxs - self.domain_mins).map(|e| e as u32)
    // }

    /// The size of each voxel part this [`Voxels`] shape.
    pub fn voxel_size(&self) -> Vector<Real> {
        self.voxel_size
    }

    /// Scale this shape.
    pub fn scaled(mut self, scale: &Vector<Real>) -> Self {
        self.voxel_size.component_mul_assign(scale);
        self
    }

    /// A reference to the chunk with id `chunk_id`.
    ///
    /// Panics if the chunk doesn’t exist.
    pub fn chunk_ref(&self, chunk_id: u32) -> VoxelsChunkRef<'_> {
        VoxelsChunkRef {
            my_id: chunk_id as usize,
            parent: self,
            states: &self.chunks[chunk_id as usize].states,
            key: &self.chunk_keys[chunk_id as usize],
        }
    }

    /// The AABB of the voxel with the given quantized `key`.
    pub fn voxel_aabb(&self, key: Point<i32>) -> Aabb {
        let center = self.voxel_center(key);
        let hext = self.voxel_size / 2.0;
        Aabb::from_half_extents(center, hext)
    }

    /// Returns the state of a given voxel.
    pub fn voxel_state(&self, key: Point<i32>) -> Option<VoxelState> {
        let vid = self.linear_index(key)?;
        Some(self.chunks[vid.chunk_id].states[vid.id_in_chunk])
    }

    /// Calculates the grid coordinates of the voxel containing the given `point`, regardless
    /// of whether this voxel is filled oor empty.
    pub fn voxel_at_point(&self, point: Point<Real>) -> Point<i32> {
        point
            .coords
            .component_div(&self.voxel_size)
            .map(|x| x.floor() as i32)
            .into()
    }

    /// Gets the voxel at the given flat voxel index.
    pub fn voxel_at_flat_id(&self, id: u32) -> Option<Point<i32>> {
        let vid = VoxelIndex::from_flat_id(id as usize);
        let chunk_key = self.chunk_keys.get(vid.chunk_id)?;
        if *chunk_key == VoxelsChunk::INVALID_CHUNK_KEY {
            return None;
        }

        Some(VoxelsChunk::voxel_key_at_id(
            *chunk_key,
            vid.id_in_chunk as u32,
        ))
    }

    /// The range of grid coordinates of voxels intersecting the given AABB.
    ///
    /// The returned range covers both empty and non-empty voxels, and is not limited to the
    /// bounds defined by [`Self::domain`].
    /// The range is semi, open, i.e., the range along each dimension `i` is understood as
    /// the semi-open interval: `range[0][i]..range[1][i]`.
    pub fn voxel_range_intersecting_local_aabb(&self, aabb: &Aabb) -> [Point<i32>; 2] {
        let mins = aabb
            .mins
            .coords
            .component_div(&self.voxel_size)
            .map(|x| x.floor() as i32);
        let maxs = aabb
            .maxs
            .coords
            .component_div(&self.voxel_size)
            .map(|x| x.ceil() as i32);
        [mins.into(), maxs.into()]
    }

    /// The AABB of a given range of voxels.
    ///
    /// The AABB is computed independently of [`Self::domain`] and independently of whether
    /// the voxels contained within are empty or not.
    pub fn voxel_range_aabb(&self, mins: Point<i32>, maxs: Point<i32>) -> Aabb {
        Aabb {
            mins: mins
                .cast::<Real>()
                .coords
                .component_mul(&self.voxel_size)
                .into(),
            maxs: maxs
                .cast::<Real>()
                .coords
                .component_mul(&self.voxel_size)
                .into(),
        }
    }

    /// Aligns the given AABB with the voxelized grid.
    ///
    /// The aligned is calculated such that the returned AABB has corners lying at the grid
    /// intersections (i.e. matches voxel corners) and fully contains the input `aabb`.
    pub fn align_aabb_to_grid(&self, aabb: &Aabb) -> Aabb {
        let mins = aabb
            .mins
            .coords
            .zip_map(&self.voxel_size, |m, sz| (m / sz).floor() * m)
            .into();
        let maxs = aabb
            .maxs
            .coords
            .zip_map(&self.voxel_size, |m, sz| (m / sz).ceil() * m)
            .into();
        Aabb { mins, maxs }
    }

    /// Iterates through every voxel intersecting the given aabb.
    ///
    /// Returns the voxel’s linearized id, center, and state.
    pub fn voxels_intersecting_local_aabb(
        &self,
        aabb: &Aabb,
    ) -> impl Iterator<Item = VoxelData> + '_ {
        let [mins, maxs] = self.voxel_range_intersecting_local_aabb(aabb);
        self.voxels_in_range(mins, maxs)
    }

    /// The center point of all the voxels in this shape (including empty ones).
    ///
    /// The voxel data associated to each center is provided to determine what kind of voxel
    /// it is (and, in particular, if it is empty or full).
    pub fn voxels(&self) -> impl Iterator<Item = VoxelData> + '_ {
        let aabb = self.chunk_bvh.root_aabb();
        self.voxels_in_range(
            self.voxel_at_point(aabb.mins),
            self.voxel_at_point(aabb.maxs),
        )
    }

    /// Iterate through the data of all the voxels within the given (semi-open) voxel grid indices.
    ///
    /// Note that this yields both empty and non-empty voxels within the range. This does not
    /// include any voxel that falls outside [`Self::domain`].
    pub fn voxels_in_range(
        &self,
        mins: Point<i32>,
        maxs: Point<i32>,
    ) -> impl Iterator<Item = VoxelData> + '_ {
        let range_aabb = Aabb::new(self.voxel_center(mins), self.voxel_center(maxs));

        self.chunk_bvh
            .leaves(move |node: &BvhNode| node.aabb().intersects(&range_aabb))
            .flat_map(move |chunk_id| {
                let chunk = self.chunk_ref(chunk_id);
                chunk.voxels_in_range(mins, maxs)
            })
    }

    fn voxel_to_chunk_key(voxel_key: Point<i32>) -> Point<i32> {
        fn div_floor(a: i32, b: usize) -> i32 {
            let sign = (a < 0) as i32;
            (a + sign) / b as i32 - sign
        }

        voxel_key.map(|e| div_floor(e, VoxelsChunk::VOXELS_PER_CHUNK_DIM))
    }

    /// Given a voxel key, returns the key of the voxel chunk that contains it, as well as the
    /// linear index of the voxel within that chunk.
    #[cfg(feature = "dim2")]
    pub(super) fn chunk_key_and_id_in_chunk(voxel_key: Point<i32>) -> (Point<i32>, usize) {
        let chunk_key = Self::voxel_to_chunk_key(voxel_key);
        // NOTE: always positive since we subtracted the smallest possible key on that chunk.
        let voxel_key_in_chunk = voxel_key - chunk_key * VoxelsChunk::VOXELS_PER_CHUNK_DIM as i32;
        let id_in_chunk = (voxel_key_in_chunk.x
            + voxel_key_in_chunk.y * VoxelsChunk::VOXELS_PER_CHUNK_DIM as i32)
            as usize;
        (chunk_key, id_in_chunk)
    }

    /// Given a voxel key, returns the key of the voxel chunk that contains it, as well as the
    /// linear index of the voxel within that chunk.
    #[cfg(feature = "dim3")]
    pub(super) fn chunk_key_and_id_in_chunk(voxel_key: Point<i32>) -> (Point<i32>, usize) {
        let chunk_key = Self::voxel_to_chunk_key(voxel_key);
        // NOTE: always positive since we subtracted the smallest possible key on that chunk.
        let voxel_key_in_chunk = voxel_key - chunk_key * VoxelsChunk::VOXELS_PER_CHUNK_DIM as i32;
        let id_in_chunk = (voxel_key_in_chunk.x
            + voxel_key_in_chunk.y * VoxelsChunk::VOXELS_PER_CHUNK_DIM as i32
            + voxel_key_in_chunk.z
                * VoxelsChunk::VOXELS_PER_CHUNK_DIM as i32
                * VoxelsChunk::VOXELS_PER_CHUNK_DIM as i32) as usize;
        (chunk_key, id_in_chunk)
    }

    /// The linearized index associated to the given voxel key.
    pub fn linear_index(&self, voxel_key: Point<i32>) -> Option<VoxelIndex> {
        let (chunk_key, id_in_chunk) = Self::chunk_key_and_id_in_chunk(voxel_key);
        let chunk_id = self.chunk_headers.get(&chunk_key)?.id;
        Some(VoxelIndex {
            chunk_id,
            id_in_chunk,
        })
    }

    /// The center of the voxel with the given key.
    pub fn voxel_center(&self, key: Point<i32>) -> Point<Real> {
        (key.cast::<Real>() + Vector::repeat(0.5))
            .coords
            .component_mul(&self.voxel_size)
            .into()
    }
}
