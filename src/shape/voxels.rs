use crate::bounding_volume::Aabb;
use crate::math::{Point, Real, Vector, DIM};
use alloc::{vec, vec::Vec};

/// The primitive shape all voxels from a [`Voxels`] is given.
#[derive(Copy, Clone, Debug, PartialEq, Default)]
#[cfg_attr(feature = "serde-serialize", derive(Serialize, Deserialize))]
pub enum VoxelPrimitiveGeometry {
    /// Each voxel is modeled as a pseudo-ball, i.e., in flat areas it will act like a planar
    /// surface but corners and edges will be rounded like a sphere.
    ///
    /// This is an approximation that is particularly relevant if the voxels are small and in large
    /// number. This can significantly improve collision-detection performances (as well as solver
    /// performances by generating less contacts points). However,this can introduce visual
    /// artifacts around edges and corners where objects in contact with the voxel will appear to
    /// slightly penetrate the corners/edges due to the spherical approximation.
    PseudoBall,
    /// Each voxel is modeled as a pseudo-cube, i.e., in flat areas it will act like a planar
    /// surface but corner and edges will be sharp like a cube.
    ///
    /// This is what you would expect for the collision to match the rendered voxels. Use
    /// [`VoxelPrimitiveGeometry::PseudoBall`] instead if some approximation around corners are acceptable
    /// in exchange for better performances.
    #[default]
    PseudoCube,
}

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
pub struct VoxelState(u8);

impl VoxelState {
    /// The value of empty voxels.
    pub const EMPTY: VoxelState = VoxelState(EMPTY_FACE_MASK);
    /// The value of a voxel with non-empty neighbors in all directions.
    pub const INTERIOR: VoxelState = VoxelState(INTERIOR_FACE_MASK);

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
    /// [`Voxels::resize_domain`].
    pub linear_id: u32,
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
/// artifacts when another objects rolls or slides against a flat voxelized surface.
///
/// The internal storage is compact (but not sparse at the moment), storing only one byte per voxel
/// in the allowed domain. This has a generally smaller memory footprint than a mesh representation
/// of the voxels.
#[derive(Clone, Debug)]
#[cfg_attr(feature = "serde-serialize", derive(Serialize, Deserialize))]
pub struct Voxels {
    domain_mins: Point<i32>,
    domain_maxs: Point<i32>,
    states: Vec<VoxelState>, // Somehow switch to a sparse representation?
    primitive_geometry: VoxelPrimitiveGeometry,
    voxel_size: Vector<Real>,
}

impl Voxels {
    /// Initializes a voxel shapes from the voxels grid coordinates.
    ///
    /// Each voxel will have its bottom-left-back corner located at
    /// `grid_coordinates * voxel_size`; and its center at `(grid_coordinates + 0.5) * voxel_size`.
    pub fn new(
        primitive_geometry: VoxelPrimitiveGeometry,
        voxel_size: Vector<Real>,
        grid_coordinates: &[Point<i32>],
    ) -> Self {
        // Ensure pseudo-balls always use uniform voxel sizes.
        let voxel_size = match primitive_geometry {
            VoxelPrimitiveGeometry::PseudoBall => Vector::repeat(voxel_size.x),
            VoxelPrimitiveGeometry::PseudoCube => voxel_size,
        };

        let mut domain_mins = grid_coordinates[0];
        let mut domain_maxs = grid_coordinates[0];

        for vox in grid_coordinates {
            domain_mins = domain_mins.inf(vox);
            domain_maxs = domain_maxs.sup(vox);
        }

        domain_maxs += Vector::repeat(1);
        let dimensions = domain_maxs - domain_mins;
        let voxels_count = dimensions.product();
        let mut result = Self {
            domain_mins,
            domain_maxs,
            states: vec![VoxelState::EMPTY; voxels_count as usize],
            primitive_geometry,
            voxel_size,
        };

        for vox in grid_coordinates {
            let index = result.linear_index(*vox);
            result.states[index as usize] = VoxelState::INTERIOR;
        }

        result.recompute_voxels_data();
        result
    }

    /// Computes a voxels shape from the set of `points`.
    ///
    /// The points are mapped to a regular grid centered at the provided point with smallest
    /// coordinates, and with grid cell size equal to `scale`. It is OK if multiple points
    /// fall into the same grid cell.
    pub fn from_points(
        primitive_geometry: VoxelPrimitiveGeometry,
        voxel_size: Vector<Real>,
        points: &[Point<Real>],
    ) -> Self {
        // Ensure pseudo-balls always use uniform voxel sizes.
        let voxel_size = match primitive_geometry {
            VoxelPrimitiveGeometry::PseudoBall => Vector::repeat(voxel_size.x),
            VoxelPrimitiveGeometry::PseudoCube => voxel_size,
        };

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
        Self::new(primitive_geometry, voxel_size, &voxels)
    }

    // TODO: support a crate like get_size2 (will require support on nalgebra too)?
    /// An approximation of the memory usage (in bytes) for this struct plus
    /// the memory it allocates dynamically.
    pub fn total_memory_size(&self) -> usize {
        size_of::<Self>() + self.heap_memory_size()
    }

    /// An approximation of the memory dynamically-allocated by this struct.
    pub fn heap_memory_size(&self) -> usize {
        // NOTE: if a new field is added to `Self`, adjust this function result.
        let Self {
            domain_mins: _,
            domain_maxs: _,
            states: data,
            primitive_geometry: _,
            voxel_size: _,
        } = self;
        data.capacity() * size_of::<VoxelState>()
    }

    /// The extents of the total axis-aligned volume covered by this [`Voxels`] shape.
    ///
    /// This accounts for all the voxels reserved in the internal buffer of `self`, including empty
    /// ones.
    pub fn extents(&self) -> Vector<Real> {
        self.dimensions()
            .cast::<Real>()
            .component_mul(&self.voxel_size)
    }

    /// The center of this shape’s domain (accounting for both empty and filled voxels).
    pub fn domain_center(&self) -> Point<Real> {
        (self
            .domain_mins
            .coords
            .cast::<Real>()
            .component_mul(&self.voxel_size)
            + self.extents() / 2.0)
            .into()
    }

    /// Sets the size of each voxel along each local coordinate axis.
    ///
    /// If [`Self::primitive_geometry`] is [`VoxelPrimitiveGeometry::PseudoBall`], then all voxels
    /// must be square, and only `size.x` is taken into account for setting the size.
    pub fn set_voxel_size(&mut self, size: Vector<Real>) {
        match self.primitive_geometry {
            VoxelPrimitiveGeometry::PseudoBall => {
                self.voxel_size = Vector::repeat(size.x);
            }
            VoxelPrimitiveGeometry::PseudoCube => {
                self.voxel_size = size;
            }
        }
    }

    /// The valid semi-open range of voxel grid indices.
    ///
    /// With `let [mins, maxs] = voxels.domain();` the valid indices along the dimension `i` are
    /// all the indices in the range `mins[i]..maxs[i]` (i.e. `maxs[i]` is excluded).
    pub fn domain(&self) -> [&Point<i32>; 2] {
        [&self.domain_mins, &self.domain_maxs]
    }

    /// The number of voxels along each coordinate axis.
    pub fn dimensions(&self) -> Vector<u32> {
        (self.domain_maxs - self.domain_mins).map(|e| e as u32)
    }

    /// The size of each voxel part this [`Voxels`] shape.
    pub fn voxel_size(&self) -> Vector<Real> {
        self.voxel_size
    }

    /// The shape each voxel is assumed to have.
    pub fn primitive_geometry(&self) -> VoxelPrimitiveGeometry {
        self.primitive_geometry
    }

    fn recompute_voxels_data(&mut self) {
        for i in 0..self.states.len() {
            let key = self.voxel_at_id(i as u32);
            self.states[i] = self.compute_voxel_state(key);
        }
    }

    /// Scale this shape.
    pub fn scaled(mut self, scale: &Vector<Real>) -> Option<Self> {
        self.voxel_size.component_mul_assign(scale);
        Some(self)
    }

    /// Sets the voxel at the given grid coordinates, returning `None` if it lies outside [`Self::domain`].
    ///
    /// See [`Self::set_voxel`] for a method that automatically resizes the internal
    /// storage of `self` if the key is out of the valid bounds.
    pub fn try_set_voxel(&mut self, key: Point<i32>, is_filled: bool) -> Option<VoxelState> {
        if key[0] < self.domain_mins[0]
            || key[0] >= self.domain_maxs[0]
            || key[1] < self.domain_mins[1]
            || key[1] >= self.domain_maxs[1]
        {
            return None;
        }

        #[cfg(feature = "dim3")]
        if key[2] < self.domain_mins[2] || key[2] >= self.domain_maxs[2] {
            return None;
        }

        let id = self.linear_index(key);
        let prev = self.states[id as usize];
        let new = if is_filled {
            VoxelState::INTERIOR
        } else {
            VoxelState::EMPTY
        };

        if prev.is_empty() ^ new.is_empty() {
            self.states[id as usize] = new;
            self.update_voxel_and_neighbors_state(key);
        }

        Some(prev)
    }

    /// Inserts a voxel at the given key, even if it is out of the bounds of this shape.
    ///
    /// If `is_filed` is `true` and the key lies out of the bounds on this shape, the internal
    /// voxels storage will be resized automatically. If a resize happens, the cost of the insertion
    /// is `O(n)` where `n` is the capacity of `self`. If no resize happens, then the cost of
    /// insertion is `O(1)`.
    ///
    /// Use [`Self::try_set_voxel`] instead for a version that will be a no-op if the provided
    /// coordinates are outside the [`Self::domain`], avoiding potential internal reallocations.
    pub fn set_voxel(&mut self, key: Point<i32>, is_filled: bool) -> Option<VoxelState> {
        if !self.is_voxel_in_bounds(key) && is_filled {
            let dims = self.dimensions();

            // Add 10% extra padding.
            let extra = dims.map(|k| k * 10 / 100);
            let mut new_domain_mins = self.domain_mins;
            let mut new_domain_maxs = self.domain_maxs;

            for k in 0..DIM {
                if key[k] < self.domain_mins[k] {
                    new_domain_mins[k] = key[k] - extra[k] as i32;
                }

                if key[k] >= self.domain_maxs[k] {
                    new_domain_maxs[k] = key[k] + extra[k] as i32 + 1;
                }
            }

            self.resize_domain(new_domain_mins, new_domain_maxs);

            self.set_voxel(key, is_filled)
        } else {
            self.set_voxel(key, is_filled)
        }
    }

    /// Set the model domain.
    ///
    /// The domain can be smaller or larger than the current one. Voxels in any overlap between
    /// the current and new domain will be preserved.
    ///
    /// If for any index `i`, `domain_maxs[i] <= domain_mins[i]`, then the new domain is invalid
    /// and this operation will result in a no-op.
    pub fn resize_domain(&mut self, domain_mins: Point<i32>, domain_maxs: Point<i32>) {
        if self.domain_mins == domain_mins && self.domain_maxs == domain_maxs {
            // Nothing to change.
            return;
        }

        if let Some(new_shape) = self.with_resized_domain(domain_mins, domain_maxs) {
            *self = new_shape;
        }
    }

    /// Clone this voxels shape with a new size.
    ///
    /// The domain can be smaller or larger than the current one. Voxels in any overlap between
    /// the current and new domain will be preserved.
    ///
    /// If for any index `i`, `domain_maxs[i] <= domain_mins[i]`, then the new domain is invalid
    /// and this operation returns `None`.
    #[must_use]
    pub fn with_resized_domain(
        &self,
        domain_mins: Point<i32>,
        domain_maxs: Point<i32>,
    ) -> Option<Self> {
        if self.domain_mins == domain_mins && self.domain_maxs == domain_maxs {
            // Nothing to change, just clone as-is.
            return Some(self.clone());
        }

        let new_dim = domain_maxs - domain_mins;
        if new_dim.iter().any(|d| *d <= 0) {
            log::error!("Invalid domain provided for resizing a voxels shape. New domain: {:?} - {:?}; new domain size: {:?}", domain_mins, domain_maxs, new_dim);
            return None;
        }

        let new_len = new_dim.iter().map(|x| *x as usize).product();

        let mut new_shape = Self {
            domain_mins,
            domain_maxs,
            states: vec![VoxelState::EMPTY; new_len],
            primitive_geometry: self.primitive_geometry,
            voxel_size: self.voxel_size,
        };

        for i in 0..self.states.len() {
            let key = self.voxel_at_id(i as u32);
            let new_i = new_shape.linear_index(key);
            new_shape.states[new_i as usize] = self.states[i];
        }

        Some(new_shape)
    }

    /// Checks if the given key is within [`Self::domain`].
    #[cfg(feature = "dim2")]
    pub fn is_voxel_in_bounds(&self, key: Point<i32>) -> bool {
        key[0] >= self.domain_mins[1]
            && key[0] < self.domain_maxs[0]
            && key[1] >= self.domain_mins[1]
            && key[1] < self.domain_maxs[1]
    }

    /// Checks if the given key is within [`Self::domain`].
    #[cfg(feature = "dim3")]
    pub fn is_voxel_in_bounds(&self, key: Point<i32>) -> bool {
        key[0] >= self.domain_mins[0]
            && key[0] < self.domain_maxs[0]
            && key[1] >= self.domain_mins[1]
            && key[1] < self.domain_maxs[1]
            && key[2] >= self.domain_mins[2]
            && key[2] < self.domain_maxs[2]
    }

    fn update_voxel_and_neighbors_state(&mut self, key: Point<i32>) {
        let key_id = self.linear_index(key) as usize;
        let mut key_data = 0;
        let center_is_empty = self.states[key_id].is_empty();

        for k in 0..DIM {
            if key[k] > self.domain_mins[k] {
                let mut left = key;
                left[k] -= 1;
                let left_id = self.linear_index(left) as usize;

                if !self.states[left_id].is_empty() {
                    if center_is_empty {
                        self.states[left_id].0 &= !(1 << (k * 2));
                    } else {
                        self.states[left_id].0 |= 1 << (k * 2);
                        key_data |= 1 << (k * 2 + 1);
                    }
                }
            }

            if key[k] + 1 < self.domain_maxs[k] {
                let mut right = key;
                right[k] += 1;
                let right_id = self.linear_index(right) as usize;

                if !self.states[right_id].is_empty() {
                    if center_is_empty {
                        self.states[right_id].0 &= !(1 << (k * 2 + 1));
                    } else {
                        self.states[right_id].0 |= 1 << (k * 2 + 1);
                        key_data |= 1 << (k * 2);
                    }
                }
            }
        }

        if !center_is_empty {
            self.states[key_id] = VoxelState(key_data);
        }
    }

    /// The AABB of the voxel with the given quantized `key`.
    pub fn voxel_aabb(&self, key: Point<i32>) -> Aabb {
        let center = self.voxel_center(key);
        let hext = self.voxel_size / 2.0;
        Aabb::from_half_extents(center, hext)
    }

    /// Returns the state of a given voxel.
    ///
    /// Panics if the key is out of the bounds defined by [`Self::domain`].
    pub fn voxel_state(&self, key: Point<i32>) -> VoxelState {
        self.states[self.linear_index(key) as usize]
    }

    /// Calculates the grid coordinates of the voxel containing the given `point`, regardless
    /// of [`Self::domain`].
    pub fn voxel_at_point_unchecked(&self, point: Point<Real>) -> Point<i32> {
        point
            .coords
            .component_div(&self.voxel_size)
            .map(|x| x.floor() as i32)
            .into()
    }

    /// Gets the key of the voxel containing the given `pt`.
    ///
    /// Note that the returned key might address a voxel that is empty.
    /// `None` is returned if the point is out of the domain of `self`.
    pub fn voxel_at_point(&self, pt: Point<Real>) -> Option<Point<i32>> {
        let quant = self.voxel_at_point_unchecked(pt);
        if quant[0] < self.domain_mins[0]
            || quant[1] < self.domain_mins[1]
            || quant[0] >= self.domain_maxs[0]
            || quant[1] >= self.domain_maxs[1]
        {
            return None;
        }

        #[cfg(feature = "dim3")]
        if quant[2] < self.domain_mins[2] || quant[2] >= self.domain_maxs[2] {
            return None;
        }

        Some(quant)
    }

    /// Clamps an arbitrary voxel into the valid domain of `self`.
    pub fn clamp_voxel(&self, key: Point<i32>) -> Point<i32> {
        key.coords
            .zip_zip_map(
                &self.domain_mins.coords,
                &self.domain_maxs.coords,
                |k, min, max| k.clamp(min, max - 1),
            )
            .into()
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
        self.voxels_in_range(self.domain_mins, self.domain_maxs)
    }

    /// Splits this voxels shape into two subshapes.
    ///
    /// The first subshape contains all the voxels which centers are inside the `aabb`.
    /// The second subshape contains all the remaining voxels.
    pub fn split_with_box(&self, aabb: &Aabb) -> (Option<Self>, Option<Self>) {
        // TODO: optimize this?
        let mut in_box = vec![];
        let mut rest = vec![];
        for vox in self.voxels() {
            if !vox.state.is_empty() {
                if aabb.contains_local_point(&vox.center) {
                    in_box.push(vox.center);
                } else {
                    rest.push(vox.center);
                }
            }
        }

        let in_box = if !in_box.is_empty() {
            Some(Voxels::from_points(
                self.primitive_geometry,
                self.voxel_size,
                &in_box,
            ))
        } else {
            None
        };

        let rest = if !rest.is_empty() {
            Some(Voxels::from_points(
                self.primitive_geometry,
                self.voxel_size,
                &rest,
            ))
        } else {
            None
        };

        (in_box, rest)
    }

    /// Iterate through the data of all the voxels within the given (semi-open) voxel grid indices.
    ///
    /// Note that this yields both empty and non-empty voxels within the range. This does not
    /// include any voxel that falls outside [`Self::domain`].
    #[cfg(feature = "dim2")]
    pub fn voxels_in_range(
        &self,
        mins: Point<i32>,
        maxs: Point<i32>,
    ) -> impl Iterator<Item = VoxelData> + '_ {
        let mins = mins.coords.sup(&self.domain_mins.coords);
        let maxs = maxs.coords.inf(&self.domain_maxs.coords);

        (mins[0]..maxs[0]).flat_map(move |ix| {
            (mins[1]..maxs[1]).map(move |iy| {
                let grid_coords = Point::new(ix, iy);
                let vid = self.linear_index(grid_coords);
                let center =
                    Vector::new(ix as Real + 0.5, iy as Real + 0.5).component_mul(&self.voxel_size);
                VoxelData {
                    linear_id: vid,
                    grid_coords,
                    center: center.into(),
                    state: self.states[vid as usize],
                }
            })
        })
    }

    /// Iterate through the data of all the voxels within the given (semi-open) voxel grid indices.
    ///
    /// Note that this yields both empty and non-empty voxels within the range. This does not
    /// include any voxel that falls outside [`Self::domain`].
    #[cfg(feature = "dim3")]
    pub fn voxels_in_range(
        &self,
        mins: Point<i32>,
        maxs: Point<i32>,
    ) -> impl Iterator<Item = VoxelData> + '_ {
        let mins = mins.coords.sup(&self.domain_mins.coords);
        let maxs = maxs.coords.inf(&self.domain_maxs.coords);

        (mins[0]..maxs[0]).flat_map(move |ix| {
            (mins[1]..maxs[1]).flat_map(move |iy| {
                (mins[2]..maxs[2]).map(move |iz| {
                    let grid_coords = Point::new(ix, iy, iz);
                    let vid = self.linear_index(grid_coords);
                    let center = Vector::new(ix as Real + 0.5, iy as Real + 0.5, iz as Real + 0.5)
                        .component_mul(&self.voxel_size)
                        .into();
                    VoxelData {
                        linear_id: vid,
                        grid_coords,
                        center,
                        state: self.states[vid as usize],
                    }
                })
            })
        })
    }

    /// The linearized index associated to the given voxel key.
    #[cfg(feature = "dim2")]
    pub fn linear_index(&self, voxel_key: Point<i32>) -> u32 {
        let dims = self.dimensions();
        let rel_key = voxel_key - self.domain_mins;
        (rel_key.x + rel_key.y * dims[0] as i32) as u32
    }

    /// The linearized index associated to the given voxel key.
    #[cfg(feature = "dim3")]
    pub fn linear_index(&self, voxel_key: Point<i32>) -> u32 {
        let dims = self.dimensions();
        let rel_key = voxel_key - self.domain_mins;
        rel_key.x as u32 + rel_key.y as u32 * dims[0] + rel_key.z as u32 * dims[0] * dims[1]
    }

    /// The key of the voxel at the given linearized index.
    #[cfg(feature = "dim2")]
    pub fn voxel_at_id(&self, linear_index: u32) -> Point<i32> {
        let dim0 = self.domain_maxs[0] - self.domain_mins[0];
        let y = linear_index as i32 / dim0;
        let x = linear_index as i32 % dim0;
        self.domain_mins + Vector::new(x, y)
    }

    /// The key of the voxel at the given linearized index.
    #[cfg(feature = "dim3")]
    pub fn voxel_at_id(&self, linear_index: u32) -> Point<i32> {
        let dims = self.dimensions();

        let d0d1 = dims[0] * dims[1];
        let z = linear_index / d0d1;
        let y = (linear_index - z * d0d1) / dims[0];
        let x = linear_index % dims[0];
        self.domain_mins + Vector::new(x as i32, y as i32, z as i32)
    }

    /// The center of the voxel with the given key.
    pub fn voxel_center(&self, key: Point<i32>) -> Point<Real> {
        (key.cast::<Real>() + Vector::repeat(0.5))
            .coords
            .component_mul(&self.voxel_size)
            .into()
    }

    fn compute_voxel_state(&self, key: Point<i32>) -> VoxelState {
        if self.states[self.linear_index(key) as usize].is_empty() {
            return VoxelState::EMPTY;
        }

        let mut occupied_faces = 0;

        for k in 0..DIM {
            let (mut prev, mut next) = (key, key);
            prev[k] -= 1;
            next[k] += 1;

            if key[k] + 1 < self.domain_maxs[k]
                && !self.states[self.linear_index(next) as usize].is_empty()
            {
                occupied_faces |= 1 << (k * 2);
            }
            if key[k] > self.domain_mins[k]
                && !self.states[self.linear_index(prev) as usize].is_empty()
            {
                occupied_faces |= 1 << (k * 2 + 1);
            }
        }

        VoxelState(occupied_faces)
    }
}

// NOTE: this code is used to generate the constant tables
// FACES_TO_VOXEL_TYPES, FACES_TO_FEATURE_MASKS, FACES_TO_OCTANT_MASKS.
#[allow(dead_code)]
#[cfg(feature = "dim2")]
#[cfg(test)]
fn gen_const_tables() {
    // The `j-th` bit of `faces_adj_to_vtx[i]` is set to 1, if the j-th face of the AABB (based on
    // the face order depicted in `AABB::FACES_VERTEX_IDS`) is adjacent to the `i` vertex of the AABB
    // (vertices are indexed as per the diagram depicted in the `FACES_VERTEX_IDS` doc.
    // Each entry of this will always have exactly 3 bits set.
    let mut faces_adj_to_vtx = [0usize; 4];

    for fid in 0..4 {
        let vids = Aabb::FACES_VERTEX_IDS[fid];
        let key = 1 << fid;
        faces_adj_to_vtx[vids.0] |= key;
        faces_adj_to_vtx[vids.1] |= key;
    }

    /*
     * FACES_TO_VOXEL_TYPES
     */
    std::println!("const FACES_TO_VOXEL_TYPES: [VoxelType; 17] = [");
    'outer: for i in 0usize..16 {
        // If any vertex of the voxel has three faces with no adjacent voxels,
        // then the voxel type is Vertex.
        for adjs in faces_adj_to_vtx.iter() {
            if (*adjs & i) == 0 {
                std::println!("VoxelType::Vertex,");
                continue 'outer;
            }
        }

        // If one face doesn’t have any adjacent voxel,
        // then the voxel type is Face.
        for fid in 0..4 {
            if ((1 << fid) & i) == 0 {
                std::println!("VoxelType::Face,");
                continue 'outer;
            }
        }
    }

    // Add final entries for special values.
    std::println!("VoxelType::Interior,");
    std::println!("VoxelType::Empty,");
    std::println!("];");

    /*
     * FACES_TO_FEATURE_MASKS
     */
    std::println!("const FACES_TO_FEATURE_MASKS: [u16; 17] = [");
    for i in 0usize..16 {
        // Each bit set indicates a convex vertex that can lead to collisions.
        // The result will be nonzero only for `VoxelType::Vertex` voxels.
        let mut vtx_key = 0;
        for (vid, adjs) in faces_adj_to_vtx.iter().enumerate() {
            if (*adjs & i) == 0 {
                vtx_key |= 1 << vid;
            }
        }

        if vtx_key != 0 {
            std::println!("0b{:b},", vtx_key as u16);
            continue;
        }

        // Each bit set indicates an exposed face that can lead to collisions.
        // The result will be nonzero only for `VoxelType::Face` voxels.
        let mut face_key = 0;
        for fid in 0..4 {
            if ((1 << fid) & i) == 0 {
                face_key |= 1 << fid;
            }
        }

        if face_key != 0 {
            std::println!("0b{:b},", face_key as u16);
            continue;
        }
    }

    std::println!("0b{:b},", u16::MAX);
    std::println!("0,");
    std::println!("];");

    /*
     * Faces to octant masks.
     */
    std::println!("const FACES_TO_OCTANT_MASKS: [u32; 17] = [");
    for i in 0usize..16 {
        // First test if we have vertices.
        let mut octant_mask = 0;
        let mut set_mask = |mask, octant| {
            // NOTE: we don’t overwrite any mask already set for the octant.
            if (octant_mask >> (octant * 3)) & 0b0111 == 0 {
                octant_mask |= mask << (octant * 3);
            }
        };

        for (vid, adjs) in faces_adj_to_vtx.iter().enumerate() {
            if (*adjs & i) == 0 {
                set_mask(1, vid);
            }
        }

        // This is the index of the normal of the faces given by
        // Aabb::FACES_VERTEX_IDS.
        const FX: u32 = OctantPattern::FACE_X;
        const FY: u32 = OctantPattern::FACE_Y;
        const FACE_NORMALS: [u32; 4] = [FX, FX, FY, FY];

        #[allow(clippy::needless_range_loop)]
        for fid in 0..4 {
            if ((1 << fid) & i) == 0 {
                let vid = Aabb::FACES_VERTEX_IDS[fid];
                let mask = FACE_NORMALS[fid];

                set_mask(mask, vid.0);
                set_mask(mask, vid.1);
            }
        }
        std::println!("0b{:b},", octant_mask);
    }
    std::println!("0,");
    std::println!("];");
}

// NOTE: this code is used to generate the constant tables
// FACES_TO_VOXEL_TYPES, FACES_TO_FEATURE_MASKS, FACES_TO_OCTANT_MASKS.
#[allow(dead_code)]
#[cfg(feature = "dim3")]
#[cfg(test)]
fn gen_const_tables() {
    // The `j-th` bit of `faces_adj_to_vtx[i]` is set to 1, if the j-th face of the AABB (based on
    // the face order depicted in `AABB::FACES_VERTEX_IDS`) is adjacent to the `i` vertex of the AABB
    // (vertices are indexed as per the diagram depicted in the `FACES_VERTEX_IDS` doc.
    // Each entry of this will always have exactly 3 bits set.
    let mut faces_adj_to_vtx = [0usize; 8];

    // The `j-th` bit of `faces_adj_to_vtx[i]` is set to 1, if the j-th edge of the AABB (based on
    // the edge order depicted in `AABB::EDGES_VERTEX_IDS`) is adjacent to the `i` vertex of the AABB
    // (vertices are indexed as per the diagram depicted in the `FACES_VERTEX_IDS` doc.
    // Each entry of this will always have exactly 2 bits set.
    let mut faces_adj_to_edge = [0usize; 12];

    for fid in 0..6 {
        let vids = Aabb::FACES_VERTEX_IDS[fid];
        let key = 1 << fid;
        faces_adj_to_vtx[vids.0] |= key;
        faces_adj_to_vtx[vids.1] |= key;
        faces_adj_to_vtx[vids.2] |= key;
        faces_adj_to_vtx[vids.3] |= key;
    }

    #[allow(clippy::needless_range_loop)]
    for eid in 0..12 {
        let evids = Aabb::EDGES_VERTEX_IDS[eid];
        for fid in 0..6 {
            let fvids = Aabb::FACES_VERTEX_IDS[fid];
            if (fvids.0 == evids.0
                || fvids.1 == evids.0
                || fvids.2 == evids.0
                || fvids.3 == evids.0)
                && (fvids.0 == evids.1
                    || fvids.1 == evids.1
                    || fvids.2 == evids.1
                    || fvids.3 == evids.1)
            {
                let key = 1 << fid;
                faces_adj_to_edge[eid] |= key;
            }
        }
    }

    /*
     * FACES_TO_VOXEL_TYPES
     */
    std::println!("const FACES_TO_VOXEL_TYPES: [VoxelType; 65] = [");
    'outer: for i in 0usize..64 {
        // If any vertex of the voxel has three faces with no adjacent voxels,
        // then the voxel type is Vertex.
        for adjs in faces_adj_to_vtx.iter() {
            if (*adjs & i) == 0 {
                std::println!("VoxelType::Vertex,");
                continue 'outer;
            }
        }

        // If any vertex of the voxel has three faces with no adjacent voxels,
        // then the voxel type is Edge.
        for adjs in faces_adj_to_edge.iter() {
            if (*adjs & i) == 0 {
                std::println!("VoxelType::Edge,");
                continue 'outer;
            }
        }

        // If one face doesn’t have any adjacent voxel,
        // then the voxel type is Face.
        for fid in 0..6 {
            if ((1 << fid) & i) == 0 {
                std::println!("VoxelType::Face,");
                continue 'outer;
            }
        }
    }

    // Add final entries for special values.
    std::println!("VoxelType::Interior,");
    std::println!("VoxelType::Empty,");
    std::println!("];");

    /*
     * FACES_TO_FEATURE_MASKS
     */
    std::println!("const FACES_TO_FEATURE_MASKS: [u16; 65] = [");
    for i in 0usize..64 {
        // Each bit set indicates a convex vertex that can lead to collisions.
        // The result will be nonzero only for `VoxelType::Vertex` voxels.
        let mut vtx_key = 0;
        for (vid, adjs) in faces_adj_to_vtx.iter().enumerate() {
            if (*adjs & i) == 0 {
                vtx_key |= 1 << vid;
            }
        }

        if vtx_key != 0 {
            std::println!("0b{:b},", vtx_key as u16);
            continue;
        }

        // Each bit set indicates a convex edge that can lead to collisions.
        // The result will be nonzero only for `VoxelType::Edge` voxels.
        let mut edge_key = 0;
        for (eid, adjs) in faces_adj_to_edge.iter().enumerate() {
            if (*adjs & i) == 0 {
                edge_key |= 1 << eid;
            }
        }

        if edge_key != 0 {
            std::println!("0b{:b},", edge_key as u16);
            continue;
        }

        // Each bit set indicates an exposed face that can lead to collisions.
        // The result will be nonzero only for `VoxelType::Face` voxels.
        let mut face_key = 0;
        for fid in 0..6 {
            if ((1 << fid) & i) == 0 {
                face_key |= 1 << fid;
            }
        }

        if face_key != 0 {
            std::println!("0b{:b},", face_key as u16);
            continue;
        }
    }

    std::println!("0b{:b},", u16::MAX);
    std::println!("0,");
    std::println!("];");

    /*
     * Faces to octant masks.
     */
    std::println!("const FACES_TO_OCTANT_MASKS: [u32; 65] = [");
    for i in 0usize..64 {
        // First test if we have vertices.
        let mut octant_mask = 0;
        let mut set_mask = |mask, octant| {
            // NOTE: we don’t overwrite any mask already set for the octant.
            if (octant_mask >> (octant * 3)) & 0b0111 == 0 {
                octant_mask |= mask << (octant * 3);
            }
        };

        for (vid, adjs) in faces_adj_to_vtx.iter().enumerate() {
            if (*adjs & i) == 0 {
                set_mask(1, vid);
            }
        }

        // This is the index of the axis porting the edges given by
        // Aabb::EDGES_VERTEX_IDS.
        const EX: u32 = OctantPattern::EDGE_X;
        const EY: u32 = OctantPattern::EDGE_Y;
        const EZ: u32 = OctantPattern::EDGE_Z;
        const EDGE_AXIS: [u32; 12] = [EX, EY, EX, EY, EX, EY, EX, EY, EZ, EZ, EZ, EZ];
        for (eid, adjs) in faces_adj_to_edge.iter().enumerate() {
            if (*adjs & i) == 0 {
                let vid = Aabb::EDGES_VERTEX_IDS[eid];
                let mask = EDGE_AXIS[eid];

                set_mask(mask, vid.0);
                set_mask(mask, vid.1);
            }
        }

        // This is the index of the normal of the faces given by
        // Aabb::FACES_VERTEX_IDS.
        const FX: u32 = OctantPattern::FACE_X;
        const FY: u32 = OctantPattern::FACE_Y;
        const FZ: u32 = OctantPattern::FACE_Z;
        const FACE_NORMALS: [u32; 6] = [FX, FX, FY, FY, FZ, FZ];

        #[allow(clippy::needless_range_loop)]
        for fid in 0..6 {
            if ((1 << fid) & i) == 0 {
                let vid = Aabb::FACES_VERTEX_IDS[fid];
                let mask = FACE_NORMALS[fid];

                set_mask(mask, vid.0);
                set_mask(mask, vid.1);
                set_mask(mask, vid.2);
                set_mask(mask, vid.3);
            }
        }
        std::println!("0b{:b},", octant_mask);
    }
    std::println!("0,");
    std::println!("];");
}

// Index to the item of FACES_TO_VOXEL_TYPES which identifies interior voxels.
#[cfg(feature = "dim2")]
const INTERIOR_FACE_MASK: u8 = 0b0000_1111;
#[cfg(feature = "dim3")]
const INTERIOR_FACE_MASK: u8 = 0b0011_1111;
// Index to the item of FACES_TO_VOXEL_TYPES which identifies empty voxels.

#[cfg(feature = "dim2")]
const EMPTY_FACE_MASK: u8 = 0b0001_0000;
#[cfg(feature = "dim3")]
const EMPTY_FACE_MASK: u8 = 0b0100_0000;

/// The voxel type deduced from adjacency information.
///
/// See the documentation of [`VoxelType`] for additional information on what each enum variant
/// means.
///
/// In 3D there are 6 neighbor faces => 64 cases + 1 empty case.
#[cfg(feature = "dim3")]
const FACES_TO_VOXEL_TYPES: [VoxelType; 65] = [
    VoxelType::Vertex,
    VoxelType::Vertex,
    VoxelType::Vertex,
    VoxelType::Edge,
    VoxelType::Vertex,
    VoxelType::Vertex,
    VoxelType::Vertex,
    VoxelType::Edge,
    VoxelType::Vertex,
    VoxelType::Vertex,
    VoxelType::Vertex,
    VoxelType::Edge,
    VoxelType::Edge,
    VoxelType::Edge,
    VoxelType::Edge,
    VoxelType::Face,
    VoxelType::Vertex,
    VoxelType::Vertex,
    VoxelType::Vertex,
    VoxelType::Edge,
    VoxelType::Vertex,
    VoxelType::Vertex,
    VoxelType::Vertex,
    VoxelType::Edge,
    VoxelType::Vertex,
    VoxelType::Vertex,
    VoxelType::Vertex,
    VoxelType::Edge,
    VoxelType::Edge,
    VoxelType::Edge,
    VoxelType::Edge,
    VoxelType::Face,
    VoxelType::Vertex,
    VoxelType::Vertex,
    VoxelType::Vertex,
    VoxelType::Edge,
    VoxelType::Vertex,
    VoxelType::Vertex,
    VoxelType::Vertex,
    VoxelType::Edge,
    VoxelType::Vertex,
    VoxelType::Vertex,
    VoxelType::Vertex,
    VoxelType::Edge,
    VoxelType::Edge,
    VoxelType::Edge,
    VoxelType::Edge,
    VoxelType::Face,
    VoxelType::Edge,
    VoxelType::Edge,
    VoxelType::Edge,
    VoxelType::Face,
    VoxelType::Edge,
    VoxelType::Edge,
    VoxelType::Edge,
    VoxelType::Face,
    VoxelType::Edge,
    VoxelType::Edge,
    VoxelType::Edge,
    VoxelType::Face,
    VoxelType::Face,
    VoxelType::Face,
    VoxelType::Face,
    VoxelType::Interior,
    VoxelType::Empty,
];

/// Indicates the convex features of each voxel that can lead to collisions.
///
/// The interpretation of each bit differs depending on the corresponding voxel type in
/// `FACES_TO_VOXEL_TYPES`:
/// - For `VoxelType::Vertex`: the i-th bit set to `1` indicates that the i-th AABB vertex is convex
///   and might lead to collisions.
/// - For `VoxelType::Edge`: the i-th bit set to `1` indicates that the i-th edge from `Aabb::EDGES_VERTEX_IDS`
///   is convex and might lead to collisions.
/// - For `VoxelType::Face`: the i-th bit set to `1` indicates that the i-th face from `Aabb::FACES_VERTEX_IDS`
///   is exposed and might lead to collisions.
#[cfg(feature = "dim3")]
const FACES_TO_FEATURE_MASKS: [u16; 65] = [
    0b11111111,
    0b10011001,
    0b1100110,
    0b1010101,
    0b110011,
    0b10001,
    0b100010,
    0b10001,
    0b11001100,
    0b10001000,
    0b1000100,
    0b1000100,
    0b10101010,
    0b10001000,
    0b100010,
    0b110000,
    0b1111,
    0b1001,
    0b110,
    0b101,
    0b11,
    0b1,
    0b10,
    0b1,
    0b1100,
    0b1000,
    0b100,
    0b100,
    0b1010,
    0b1000,
    0b10,
    0b100000,
    0b11110000,
    0b10010000,
    0b1100000,
    0b1010000,
    0b110000,
    0b10000,
    0b100000,
    0b10000,
    0b11000000,
    0b10000000,
    0b1000000,
    0b1000000,
    0b10100000,
    0b10000000,
    0b100000,
    0b10000,
    0b111100000000,
    0b100100000000,
    0b11000000000,
    0b1100,
    0b1100000000,
    0b100000000,
    0b1000000000,
    0b1000,
    0b110000000000,
    0b100000000000,
    0b10000000000,
    0b100,
    0b11,
    0b10,
    0b1,
    0b1111111111111111,
    0,
];

/// Each octant is assigned three contiguous bits.
#[cfg(feature = "dim3")]
const FACES_TO_OCTANT_MASKS: [u32; 65] = [
    0b1001001001001001001001,
    0b1010010001001010010001,
    0b10001001010010001001010,
    0b10010010010010010010010,
    0b11011001001011011001001,
    0b11111010001011111010001,
    0b111011001010111011001010,
    0b111111010010111111010010,
    0b1001011011001001011011,
    0b1010111011001010111011,
    0b10001011111010001011111,
    0b10010111111010010111111,
    0b11011011011011011011011,
    0b11111111011011111111011,
    0b111011011111111011011111,
    0b111111111111111111111111,
    0b100100100100001001001001,
    0b100110110100001010010001,
    0b110100100110010001001010,
    0b110110110110010010010010,
    0b101101100100011011001001,
    0b101000110100011111010001,
    0b101100110111011001010,
    0b110110111111010010,
    0b100100101101001001011011,
    0b100110000101001010111011,
    0b110100101000010001011111,
    0b110110000000010010111111,
    0b101101101101011011011011,
    0b101000000101011111111011,
    0b101101000111011011111,
    0b111111111111,
    0b1001001001100100100100,
    0b1010010001100110110100,
    0b10001001010110100100110,
    0b10010010010110110110110,
    0b11011001001101101100100,
    0b11111010001101000110100,
    0b111011001010000101100110,
    0b111111010010000000110110,
    0b1001011011100100101101,
    0b1010111011100110000101,
    0b10001011111110100101000,
    0b10010111111110110000000,
    0b11011011011101101101101,
    0b11111111011101000000101,
    0b111011011111000101101000,
    0b111111111111000000000000,
    0b100100100100100100100100,
    0b100110110100100110110100,
    0b110100100110110100100110,
    0b110110110110110110110110,
    0b101101100100101101100100,
    0b101000110100101000110100,
    0b101100110000101100110,
    0b110110000000110110,
    0b100100101101100100101101,
    0b100110000101100110000101,
    0b110100101000110100101000,
    0b110110000000110110000000,
    0b101101101101101101101101,
    0b101000000101101000000101,
    0b101101000000101101000,
    0b0,
    0,
];

#[cfg(feature = "dim2")]
const FACES_TO_VOXEL_TYPES: [VoxelType; 17] = [
    VoxelType::Vertex,
    VoxelType::Vertex,
    VoxelType::Vertex,
    VoxelType::Face,
    VoxelType::Vertex,
    VoxelType::Vertex,
    VoxelType::Vertex,
    VoxelType::Face,
    VoxelType::Vertex,
    VoxelType::Vertex,
    VoxelType::Vertex,
    VoxelType::Face,
    VoxelType::Face,
    VoxelType::Face,
    VoxelType::Face,
    VoxelType::Interior,
    VoxelType::Empty,
];

#[cfg(feature = "dim2")]
const FACES_TO_FEATURE_MASKS: [u16; 17] = [
    0b1111,
    0b1001,
    0b110,
    0b1100,
    0b11,
    0b1,
    0b10,
    0b1000,
    0b1100,
    0b1000,
    0b100,
    0b100,
    0b11,
    0b10,
    0b1,
    0b1111111111111111,
    0,
];

// NOTE: in 2D we are also using 3 bits per octant even though we technically only need two.
//       This keeps some collision-detection easier by avoiding some special-casing.
#[cfg(feature = "dim2")]
const FACES_TO_OCTANT_MASKS: [u32; 17] = [
    0b1001001001,
    0b1011011001,
    0b11001001011,
    0b11011011011,
    0b10010001001,
    0b10000011001,
    0b10001011,
    0b11011,
    0b1001010010,
    0b1011000010,
    0b11001010000,
    0b11011000000,
    0b10010010010,
    0b10000000010,
    0b10010000,
    0b0,
    0,
];

#[cfg(test)]
mod test {
    #[test]
    fn gen_const_tables() {
        super::gen_const_tables();
    }
}
