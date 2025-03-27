use crate::bounding_volume::Aabb;
use crate::math::{Point, Real, Vector, DIM};

/// The primitive shape all voxels from a [`Voxels`] is given.
#[derive(Copy, Clone, Debug, Default, PartialEq, Eq)]
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

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct OctantPattern;

// NOTE: it is important that the max value of any OctantPattern variant
//       is 7 because we don’t allocate more than 3 bits to store it in
//      `FACES_TO_OCTANT_MASKS`.
/// Indicates the local shape of a voxel on each octant.
///
/// This provides geometric information of the shape’s exposed features on each octant.
/// This is an alternative to `FACES_TO_FEATURE_MASKS` that can be more convenient for some
/// collision-detection algorithms.
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
pub struct VoxelData(u8);

impl VoxelData {
    /// The value of empty voxels.
    pub const EMPTY: VoxelData = VoxelData(EMPTY_FACE_MASK);
    /// The value of a voxel with non-empty neighbors in all directions.
    pub const INTERIOR: VoxelData = VoxelData(INTERIOR_FACE_MASK);

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
    pub const fn voxel_type(self) -> VoxelType {
        FACES_TO_VOXEL_TYPES[self.0 as usize]
    }

    // Bitmask indicating what vertices, edges, or faces of the voxel are "free".
    pub const fn feature_mask(self) -> u16 {
        FACES_TO_FEATURE_MASKS[self.0 as usize]
    }
    pub const fn octant_mask(self) -> u32 {
        FACES_TO_OCTANT_MASKS[self.0 as usize]
    }
}

#[derive(Clone, Debug)]
#[cfg_attr(feature = "serde-serialize", derive(Serialize, Deserialize))]
pub struct Voxels {
    pub(crate) dimensions: [u32; DIM],
    pub(crate) data: Vec<VoxelData>, // TODO: use a Vob?
    primitive_geometry: VoxelPrimitiveGeometry,
    pub scale: Real,
    pub origin: Point<Real>,
}

impl Voxels {
    /// Computes a voxels shape from the set of `points`.
    ///
    /// The points are mapped to a regular grid centered at the provided point with smallest
    /// coordinates, and with grid cell size equal to `scale`. It is OK if multiple points
    /// fall into the same grid cell.
    pub fn from_points(
        primitive_geometry: VoxelPrimitiveGeometry,
        points: &[Point<Real>],
        scale: Real,
    ) -> Self {
        let mut aabb = Aabb::from_points(points);
        aabb.mins -= Vector::repeat(scale / 2.0);
        aabb.maxs += Vector::repeat(scale / 2.0); // NOTE: is this + necessary?

        // Discretize the points on the grid.
        let dimensions = ((aabb.maxs - aabb.mins) / scale).map(|x| x.ceil() as u32);
        let len = dimensions.product();
        let mut result = Self {
            dimensions: dimensions.into(),
            data: vec![VoxelData::EMPTY; len as usize],
            primitive_geometry,
            scale,
            origin: aabb.mins,
        };

        for pt in points {
            let coords = (pt - aabb.mins).map(|x| (x / scale).floor() as u32);
            let index = result.linear_index(coords.into());
            // The precise voxel data will be computed in `recompute_voxels_data` below.
            result.data[index as usize] = VoxelData::INTERIOR;
        }

        result.recompute_voxels_data();
        result
    }

    /// The total axis-aligned volume covered by this [`Voxels`] shape.
    ///
    /// This accounts for all the voxels in `self`, including empty ones.
    pub fn extents(&self) -> Vector<Real> {
        Vector::from(self.dimensions).cast::<Real>() * self.scale
    }

    /// The size of each voxel part this [`Voxels`] shape.
    pub fn voxel_size(&self) -> Real {
        self.scale
    }

    pub fn primitive_geometry(&self) -> VoxelPrimitiveGeometry {
        self.primitive_geometry
    }

    fn recompute_voxels_data(&mut self) {
        for i in 0..self.data.len() {
            let key = self.to_key(i as u32);
            self.data[i] = self.compute_voxel_data(key);
        }
    }

    /// Scale this shape. Returns `None` if the scale is non-uniform (not supported yet).
    pub fn scaled(mut self, scale: &Vector<Real>) -> Option<Self> {
        #[cfg(feature = "dim2")]
        let scale_is_uniform = scale.x == scale.y;
        #[cfg(feature = "dim3")]
        let scale_is_uniform = scale.x == scale.y && scale.y == scale.z;
        if !scale_is_uniform {
            // TODO: what about non-uniform scale?
            None
        } else {
            self.scale *= scale.x;
            self.origin *= scale.x;
            Some(self)
        }
    }

    fn quantify_point(&self, pt: Point<Real>) -> Vector<u32> {
        ((pt - self.origin) / self.scale).map(|x| x.floor().max(0.0) as u32)
    }

    pub fn voxels_intersecting_local_aabb(
        &self,
        aabb: &Aabb,
    ) -> impl Iterator<Item = (Point<Real>, VoxelData)> + '_ {
        let dims = Vector::from(self.dimensions);
        let mins = ((aabb.mins - self.origin) / self.scale)
            .map(|x| x.floor().max(0.0) as u32)
            .inf(&dims);
        let maxs = ((aabb.maxs - self.origin) / self.scale)
            .map(|x| x.ceil().max(0.0) as u32)
            .inf(&dims);

        self.centers_range(mins.into(), maxs.into())
    }

    /// The center point of all the voxels in this shape (including empty ones).
    ///
    /// The voxel data associated to each center is provided to determine what kind of voxel
    /// it is (and, in particular, if it is empty or full).
    pub fn centers(&self) -> impl Iterator<Item = (Point<Real>, VoxelData)> + '_ {
        self.centers_range([0; DIM], self.dimensions)
    }

    /// Splits this voxels shape into two subshapes.
    ///
    /// The first subshape contains all the voxels which centers are inside the `aabb`.
    /// The second subshape contains all the remaining voxels.
    pub fn split_with_box(&self, aabb: &Aabb) -> (Option<Self>, Option<Self>) {
        // TODO: optimize this?
        let mut in_box = vec![];
        let mut rest = vec![];
        for (center, voxel) in self.centers() {
            if !voxel.is_empty() {
                if aabb.contains_local_point(&center) {
                    in_box.push(center);
                } else {
                    rest.push(center);
                }
            }
        }

        let in_box = if !in_box.is_empty() {
            Some(Voxels::from_points(
                self.primitive_geometry,
                &in_box,
                self.scale,
            ))
        } else {
            None
        };

        let rest = if !rest.is_empty() {
            Some(Voxels::from_points(
                self.primitive_geometry,
                &rest,
                self.scale,
            ))
        } else {
            None
        };

        (in_box, rest)
    }

    #[cfg(feature = "dim2")]
    fn centers_range(
        &self,
        mins: [u32; DIM],
        maxs: [u32; DIM],
    ) -> impl Iterator<Item = (Point<Real>, VoxelData)> + '_ {
        (mins[0]..maxs[0]).flat_map(move |ix| {
            (mins[1]..maxs[1]).map(move |iy| {
                let vid = self.linear_index([ix, iy]);
                let center =
                    self.origin + Vector::new(ix as Real + 0.5, iy as Real + 0.5) * self.scale;
                (center, self.data[vid as usize])
            })
        })
    }

    #[cfg(feature = "dim3")]
    fn centers_range(
        &self,
        mins: [u32; DIM],
        maxs: [u32; DIM],
    ) -> impl Iterator<Item = (Point<Real>, VoxelData)> + '_ {
        (mins[0]..maxs[0]).flat_map(move |ix| {
            (mins[1]..maxs[1]).flat_map(move |iy| {
                (mins[2]..maxs[2]).map(move |iz| {
                    let vid = self.linear_index([ix, iy, iz]);
                    let center = self.origin
                        + Vector::new(ix as Real + 0.5, iy as Real + 0.5, iz as Real + 0.5)
                            * self.scale;
                    (center, self.data[vid as usize])
                })
            })
        })
    }

    #[cfg(feature = "dim2")]
    fn linear_index(&self, voxel_id: [u32; DIM]) -> u32 {
        voxel_id[0] + voxel_id[1] * self.dimensions[0]
    }

    #[cfg(feature = "dim3")]
    fn linear_index(&self, voxel_id: [u32; DIM]) -> u32 {
        voxel_id[0]
            + voxel_id[1] * self.dimensions[0]
            + voxel_id[2] * self.dimensions[0] * self.dimensions[1]
    }

    #[cfg(feature = "dim2")]
    fn to_key(&self, linear_index: u32) -> [u32; DIM] {
        let y = linear_index / self.dimensions[0];
        let x = linear_index % self.dimensions[0];
        [x, y]
    }

    #[cfg(feature = "dim3")]
    fn to_key(&self, linear_index: u32) -> [u32; DIM] {
        let d0d1 = self.dimensions[0] * self.dimensions[1];
        let z = linear_index / d0d1;
        let y = (linear_index - z * d0d1) / self.dimensions[0];
        let x = linear_index % self.dimensions[0];
        [x, y, z]
    }

    fn compute_voxel_data(&self, key: [u32; DIM]) -> VoxelData {
        if self.data[self.linear_index(key) as usize].is_empty() {
            return VoxelData::EMPTY;
        }

        let mut occupied_faces = 0;

        for k in 0..DIM {
            let (mut prev, mut next) = (key, key);
            prev[k] = prev[k].saturating_sub(1);
            next[k] += 1;

            if key[k] < self.dimensions[k] - 1
                && !self.data[self.linear_index(next) as usize].is_empty()
            {
                occupied_faces |= 1 << (k * 2);
            }
            if key[k] > 0 && !self.data[self.linear_index(prev) as usize].is_empty() {
                occupied_faces |= 1 << (k * 2 + 1);
            }
        }

        VoxelData(occupied_faces)
    }
}

// NOTE: this code is used to generate the constant tables
// FACES_TO_VOXEL_TYPES, FACES_TO_FEATURE_MASKS, FACES_TO_OCTANT_MASKS.
#[allow(dead_code)]
#[cfg(feature = "dim2")]
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
    println!("const FACES_TO_VOXEL_TYPES: [VoxelType; 17] = [");
    'outer: for i in 0usize..16 {
        // If any vertex of the voxel has three faces with no adjacent voxels,
        // then the voxel type is Vertex.
        for adjs in faces_adj_to_vtx.iter() {
            if (*adjs & i) == 0 {
                println!("VoxelType::Vertex,");
                continue 'outer;
            }
        }

        // If one face doesn’t have any adjacent voxel,
        // then the voxel type is Face.
        for fid in 0..4 {
            if ((1 << fid) & i) == 0 {
                println!("VoxelType::Face,");
                continue 'outer;
            }
        }
    }

    // Add final entries for special values.
    println!("VoxelType::Interior,");
    println!("VoxelType::Empty,");
    println!("];");

    /*
     * FACES_TO_FEATURE_MASKS
     */
    println!("const FACES_TO_FEATURE_MASKS: [u16; 17] = [");
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
            println!("0b{:b},", vtx_key as u16);
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
            println!("0b{:b},", face_key as u16);
            continue;
        }
    }

    println!("0b{:b},", u16::MAX);
    println!("0,");
    println!("];");

    /*
     * Faces to octant masks.
     */
    println!("const FACES_TO_OCTANT_MASKS: [u32; 17] = [");
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
        for fid in 0..4 {
            if ((1 << fid) & i) == 0 {
                let vid = Aabb::FACES_VERTEX_IDS[fid];
                let mask = FACE_NORMALS[fid];

                set_mask(mask, vid.0);
                set_mask(mask, vid.1);
            }
        }
        println!("0b{:b},", octant_mask);
    }
    println!("0,");
    println!("];");
}

// NOTE: this code is used to generate the constant tables
// FACES_TO_VOXEL_TYPES, FACES_TO_FEATURE_MASKS, FACES_TO_OCTANT_MASKS.
#[allow(dead_code)]
#[cfg(feature = "dim3")]
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
    println!("const FACES_TO_VOXEL_TYPES: [VoxelType; 65] = [");
    'outer: for i in 0usize..64 {
        // If any vertex of the voxel has three faces with no adjacent voxels,
        // then the voxel type is Vertex.
        for adjs in faces_adj_to_vtx.iter() {
            if (*adjs & i) == 0 {
                println!("VoxelType::Vertex,");
                continue 'outer;
            }
        }

        // If any vertex of the voxel has three faces with no adjacent voxels,
        // then the voxel type is Edge.
        for adjs in faces_adj_to_edge.iter() {
            if (*adjs & i) == 0 {
                println!("VoxelType::Edge,");
                continue 'outer;
            }
        }

        // If one face doesn’t have any adjacent voxel,
        // then the voxel type is Face.
        for fid in 0..6 {
            if ((1 << fid) & i) == 0 {
                println!("VoxelType::Face,");
                continue 'outer;
            }
        }
    }

    // Add final entries for special values.
    println!("VoxelType::Interior,");
    println!("VoxelType::Empty,");
    println!("];");

    /*
     * FACES_TO_FEATURE_MASKS
     */
    println!("const FACES_TO_FEATURE_MASKS: [u16; 65] = [");
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
            println!("0b{:b},", vtx_key as u16);
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
            println!("0b{:b},", edge_key as u16);
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
            println!("0b{:b},", face_key as u16);
            continue;
        }
    }

    println!("0b{:b},", u16::MAX);
    println!("0,");
    println!("];");

    /*
     * Faces to octant masks.
     */
    println!("const FACES_TO_OCTANT_MASKS: [u32; 65] = [");
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
        println!("0b{:b},", octant_mask);
    }
    println!("0,");
    println!("];");
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
///                            and might lead to collisions.
/// - For `VoxelType::Edge`: the i-th bit set to `1` indicates that the i-th edge from `Aabb::EDGES_VERTEX_IDS`
///                          is convex and might lead to collisions.
/// - For `VoxelType::Face`: the i-th bit set to `1` indicates that the i-th face from `Aabb::FACES_VERTEX_IDS`
///                          is exposed and might lead to collisions.
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
