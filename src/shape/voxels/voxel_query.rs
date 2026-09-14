use crate::math::{ivect_to_vect, vect_to_ivect, IVector, Vector};

use crate::bounding_volume::Aabb;
use crate::shape::{VoxelData, VoxelState, VoxelType, Voxels};

/// Abstraction over the storage of a shape made of axis-aligned, uniformly sized voxels.
///
/// Parry's voxel collision-detection algorithms (contact manifolds, intersection tests,
/// linear and nonlinear shape-casting) are
/// written against this trait rather than against the concrete [`Voxels`] shape. Implementing
/// it for a custom sparse data-structure (chunked grid, octree, VDB-like tree, etc.) lets these
/// algorithms run directly on that structure without copying it into a [`Voxels`] shape,
/// typically by calling the generic query functions from a custom
/// [`QueryDispatcher`](crate::query::QueryDispatcher).
///
/// # Grid conventions
///
/// Voxels are identified by their integer grid coordinates `key`. The voxel with coordinates
/// `key` covers the world-space (well, shape-local-space) range
/// `[key * voxel_size, (key + 1) * voxel_size]`, so its center is at
/// `(key + 0.5) * voxel_size`. Grid ranges are always given as semi-open intervals
/// `[mins, maxs)`: `mins` is included, `maxs` is excluded.
///
/// # Voxel views and neighborhood states
///
/// Lookups and iterators don't yield a fixed data struct: they yield storage-defined voxel
/// *views* ([`Self::Voxel`], bounded by [`QueriedVoxel`]). A view exposes cheap per-voxel
/// data — grid coordinates, center, and the coarse [`QueriedVoxel::voxel_type`], which a
/// sparse storage can pack in two bits per stored voxel (with empty voxels simply absent).
///
/// Contact-manifold computation additionally needs to know *which* of a voxel's immediate
/// axis-aligned neighbors are filled — a [`VoxelState`] — to avoid hitting the "internal
/// edges" between adjacent voxels. It obtains this from [`QueriedVoxel::voxel_state`], and
/// only for the few voxels that are actual contact candidates, never during bulk iteration.
/// Since views can borrow from their storage, they can compute the state on demand from
/// local context (e.g. leaf-local reads in a sparse tree);
/// [`VoxelState::with_filled_neighbors`] builds the state from occupancy alone, while storages
/// like [`Voxels`] that persist the state (one byte per voxel) just hand out the stored value.
///
/// # Note for implementors
///
/// This trait is not dyn-compatible (`voxels_in_range` returns `impl Iterator`). The generic
/// query functions are monomorphized for each storage type. To plug a custom storage into a
/// physics pipeline, wrap it in a type implementing [`Shape`](crate::shape::Shape) (typically
/// with [`ShapeType::Custom`](crate::shape::ShapeType::Custom)) and dispatch to the generic
/// voxel query functions from a custom `QueryDispatcher`.
pub trait VoxelQuery {
    /// The view type this storage hands out for a single voxel.
    ///
    /// Views can borrow from the storage (e.g. hold a cursor into a sparse tree), letting
    /// [`QueriedVoxel::voxel_state`] read neighborhood information from local context
    /// instead of independent whole-storage lookups.
    type Voxel<'a>: QueriedVoxel<'a>
    where
        Self: 'a;

    /// The size of each voxel along each local coordinate axis.
    fn voxel_size(&self) -> Vector;

    /// The semi-open range `[mins, maxs)` of grid coordinates covered by this shape.
    ///
    /// This must be a conservative bound: every non-empty voxel must lie within the returned
    /// range, but the range may also cover empty voxels.
    fn domain(&self) -> [IVector; 2];

    /// Iterates through the voxels within the given semi-open grid coordinate range.
    ///
    /// Implementations must yield every non-empty voxel with grid coordinates in
    /// `[mins, maxs)` exactly once. They may additionally yield empty voxels within that
    /// range (callers filter on [`QueriedVoxel::voxel_type`]), but must never yield a
    /// voxel outside of the range.
    fn voxels_in_range(
        &self,
        mins: IVector,
        maxs: IVector,
    ) -> impl Iterator<Item = Self::Voxel<'_>>;

    /// Iterates through every voxel of this shape.
    ///
    /// This is equivalent to [`Self::voxels_in_range`] applied to the whole [`Self::domain`].
    fn voxels(&self) -> impl Iterator<Item = Self::Voxel<'_>> {
        let [mins, maxs] = self.domain();
        self.voxels_in_range(mins, maxs)
    }

    /// Iterates through every voxel intersecting the given local-space AABB.
    fn voxels_intersecting_local_aabb(&self, aabb: &Aabb) -> impl Iterator<Item = Self::Voxel<'_>> {
        let [mins, maxs] = self.voxel_range_intersecting_local_aabb(aabb);
        self.voxels_in_range(mins, maxs)
    }

    /// The grid coordinates of the voxel containing the given local-space point.
    ///
    /// The returned coordinates are valid regardless of whether the corresponding voxel
    /// is filled, empty, or outside of [`Self::domain`].
    fn voxel_at_point(&self, point: Vector) -> IVector {
        vect_to_ivect((point / self.voxel_size()).floor())
    }

    /// The local-space center of the voxel with the given grid coordinates.
    fn voxel_center(&self, key: IVector) -> Vector {
        (ivect_to_vect(key) + Vector::splat(0.5)) * self.voxel_size()
    }

    /// The local-space AABB of the voxel with the given grid coordinates.
    fn voxel_aabb(&self, key: IVector) -> Aabb {
        let center = self.voxel_center(key);
        Aabb::from_half_extents(center, self.voxel_size() / 2.0)
    }

    /// The semi-open range of grid coordinates of the voxels intersecting the given AABB.
    ///
    /// The returned range covers both empty and non-empty voxels, and is not limited to the
    /// bounds defined by [`Self::domain`].
    fn voxel_range_intersecting_local_aabb(&self, aabb: &Aabb) -> [IVector; 2] {
        let mins = vect_to_ivect((aabb.mins / self.voxel_size()).floor());
        let maxs = vect_to_ivect((aabb.maxs / self.voxel_size()).ceil());
        [mins, maxs]
    }

    /// The local-space AABB of the given semi-open range of voxel grid coordinates.
    fn voxel_range_aabb(&self, mins: IVector, maxs: IVector) -> Aabb {
        Aabb {
            mins: ivect_to_vect(mins) * self.voxel_size(),
            maxs: ivect_to_vect(maxs) * self.voxel_size(),
        }
    }

    /// Aligns the given AABB with the voxelized grid.
    ///
    /// The returned AABB has corners lying at the grid intersections (i.e. matches voxel
    /// corners) and fully contains the input `aabb`.
    fn align_aabb_to_grid(&self, aabb: &Aabb) -> Aabb {
        let mins = (aabb.mins / self.voxel_size()).floor() * self.voxel_size();
        let maxs = (aabb.maxs / self.voxel_size()).ceil() * self.voxel_size();
        Aabb { mins, maxs }
    }

    /// The local-space AABB of this voxels shape.
    fn local_aabb(&self) -> Aabb {
        let [mins, maxs] = self.domain();
        self.voxel_range_aabb(mins, maxs)
    }
}

/// A single voxel handed out by a [`VoxelQuery`] storage.
///
/// This is the item type of the [`VoxelQuery`] lookups and iterators. Storages define their
/// own implementor (see [`VoxelQuery::Voxel`]), which may borrow from the storage so that
/// [`Self::voxel_state`] can be computed lazily from local context.
///
/// The coarse [`Self::voxel_type`] must be cheap: it is read during bulk iteration. The full
/// [`Self::voxel_state`] is only requested for contact-candidate voxels; implementations
/// must keep the two consistent (`self.voxel_state().voxel_type() == self.voxel_type()`).
pub trait QueriedVoxel<'a> {
    /// The type of this voxel: empty, or how it is exposed on the shape's surface.
    fn voxel_type(&self) -> VoxelType;

    /// The neighborhood state of this voxel, indicating which of its immediate
    /// axis-aligned neighbors are filled.
    fn voxel_state(&self) -> VoxelState;

    /// A stable, storage-defined identifier of this voxel.
    ///
    /// For the [`Voxels`] shape this is the flattened form of [`Voxels::linear_index`].
    /// This identifier can be invalidated after the voxels shape is modified (e.g. by a call
    /// to [`Voxels::set_voxel`], or [`Voxels::crop`]).
    /// For stable references to voxels, always use [`Self::grid_coords`]. Only meaningful
    /// for non-empty voxels.
    fn linear_id(&self) -> u32;

    /// The voxel's integer grid coordinates.
    fn grid_coords(&self) -> IVector;

    /// The voxel's center position in the local-space of the voxels shape it is part of.
    fn center(&self) -> Vector;
}

impl QueriedVoxel<'_> for VoxelData {
    fn voxel_type(&self) -> VoxelType {
        self.state.voxel_type()
    }

    fn voxel_state(&self) -> VoxelState {
        self.state
    }

    fn linear_id(&self) -> u32 {
        self.linear_id
    }

    fn grid_coords(&self) -> IVector {
        self.grid_coords
    }

    fn center(&self) -> Vector {
        self.center
    }
}
impl VoxelQuery for Voxels {
    type Voxel<'a> = VoxelData;

    #[inline]
    fn voxel_size(&self) -> Vector {
        self.voxel_size()
    }

    #[inline]
    fn domain(&self) -> [IVector; 2] {
        self.domain()
    }

    #[inline]
    fn voxels_in_range(&self, mins: IVector, maxs: IVector) -> impl Iterator<Item = VoxelData> {
        self.voxels_in_range(mins, maxs)
    }

    #[inline]
    fn voxels(&self) -> impl Iterator<Item = VoxelData> {
        self.voxels()
    }

    #[inline]
    fn local_aabb(&self) -> Aabb {
        self.local_aabb()
    }
}
