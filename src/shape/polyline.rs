use crate::bounding_volume::Aabb;
use crate::math::{Pose, Vector};
use crate::partitioning::{Bvh, BvhBuildStrategy};
use crate::query::{PointProjection, PointQueryWithLocation};
use crate::shape::composite_shape::CompositeShape;
use crate::shape::{
    FeatureId, Segment, SegmentPointLocation, SegmentPseudoNormals, Shape, TypedCompositeShape,
};
#[cfg(feature = "alloc")]
use alloc::vec::Vec;

use crate::query::details::NormalConstraints;

#[cfg(feature = "dim2")]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Deserialize, rkyv::Serialize)
)]
#[repr(C)]
#[derive(Clone, Copy, Debug, Default, Eq, Hash, Ord, PartialEq, PartialOrd)]
/// Controls how a [`Polyline`] is loaded.
pub struct PolylineFlags(u8);

#[cfg(feature = "dim2")]
bitflags::bitflags! {
    impl PolylineFlags: u8 {
        /// If set, the polyline is treated as one-sided: a pseudo-normal is computed at every
        /// vertex and contact normals are clamped to the outward side, the *right* of each
        /// segment's direction. The solid must be wound counter-clockwise, so the outward side is
        /// on the right. This removes the spurious sideways push a body gets at a convex corner of
        /// a double-sided polyline. This one flag covers what `TriMesh` splits across
        /// `TriMeshFlags::ORIENTED` (compute pseudo-normals) and `TriMeshFlags::FIX_INTERNAL_EDGES`
        /// (use them to clamp contacts).
        const ORIENTED = 1;
    }
}

#[derive(Clone, Debug)]
#[cfg_attr(feature = "serde-serialize", derive(Serialize, Deserialize))]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Deserialize, rkyv::Serialize)
)]
/// A polyline shape formed by connected line segments.
///
/// A polyline is a sequence of line segments (edges) connecting vertices. It can be open
/// (not forming a closed loop) or closed (where the last vertex connects back to the first).
/// Polylines are commonly used for paths, boundaries, and 2D/3D curves.
///
/// # Structure
///
/// A polyline consists of:
/// - **Vertices**: Vectors in 2D or 3D space
/// - **Indices**: Pairs of vertex indices defining each segment
/// - **BVH**: Bounding Volume Hierarchy for fast spatial queries
///
/// # Properties
///
/// - **Composite shape**: Made up of multiple segments
/// - **1-dimensional**: Has length but no volume
/// - **Flexible topology**: Can be open or closed, branching or linear
/// - **Accelerated queries**: Uses BVH for efficient collision detection
///
/// # Use Cases
///
/// Polylines are ideal for:
/// - **Paths and roads**: Navigation paths, road networks
/// - **Terrain boundaries**: Cliff edges, coastlines, level boundaries
/// - **Outlines**: 2D shape outlines, contours
/// - **Wire frames**: Simplified representations of complex shapes
/// - **Motion paths**: Character movement paths, camera rails
///
/// # Example
///
/// ```rust
/// # #[cfg(all(feature = "dim3", feature = "f32"))] {
/// use parry3d::shape::Polyline;
/// use parry3d::math::Vector;
///
/// // Create a simple L-shaped polyline
/// let vertices = vec![
///     Vector::ZERO,
///     Vector::new(1.0, 0.0, 0.0),
///     Vector::new(1.0, 1.0, 0.0),
/// ];
///
/// // Indices are automatically generated to connect consecutive vertices
/// let polyline = Polyline::new(vertices, None);
///
/// // The polyline has 2 segments: (0,1) and (1,2)
/// assert_eq!(polyline.num_segments(), 2);
/// assert_eq!(polyline.vertices().len(), 3);
/// # }
/// ```
///
/// # Custom Connectivity
///
/// You can provide custom indices to create non-sequential connections:
///
/// ```rust
/// # #[cfg(all(feature = "dim2", feature = "f32"))] {
/// use parry2d::shape::Polyline;
/// use parry2d::math::Vector;
///
/// // Create a triangle polyline (closed loop)
/// let vertices = vec![
///     Vector::ZERO,
///     Vector::new(1.0, 0.0),
///     Vector::new(0.5, 1.0),
/// ];
///
/// // Manually specify edges to create a closed triangle
/// let indices = vec![
///     [0, 1],  // Bottom edge
///     [1, 2],  // Right edge
///     [2, 0],  // Left edge (closes the loop)
/// ];
///
/// let polyline = Polyline::new(vertices, Some(indices));
/// assert_eq!(polyline.num_segments(), 3);
/// # }
/// ```
pub struct Polyline {
    bvh: Bvh,
    vertices: Vec<Vector>,
    indices: Vec<[u32; 2]>,
    /// Per-vertex outward pseudo-normals, present when [`PolylineFlags::ORIENTED`] is set; contact
    /// normals are then clamped to one side so the polyline acts as a one-sided surface.
    #[cfg(feature = "dim2")]
    pseudo_normals: Option<Vec<Vector>>,
    #[cfg(feature = "dim2")]
    flags: PolylineFlags,
}

impl Polyline {
    /// Creates a new polyline from a vertex buffer and an optional index buffer.
    ///
    /// This is the main constructor for creating a polyline. If no indices are provided,
    /// the vertices will be automatically connected in sequence (vertex 0 to 1, 1 to 2, etc.).
    ///
    /// # Arguments
    ///
    /// * `vertices` - A vector of points defining the polyline vertices
    /// * `indices` - Optional vector of `[u32; 2]` pairs defining which vertices connect.
    ///   If `None`, vertices are connected sequentially.
    ///
    /// # Returns
    ///
    /// A new `Polyline` with an internal BVH for accelerated queries.
    ///
    /// # Example
    ///
    /// ```
    /// # #[cfg(all(feature = "dim3", feature = "f32"))] {
    /// use parry3d::shape::Polyline;
    /// use parry3d::math::Vector;
    ///
    /// // Create a zigzag path with automatic sequential connections
    /// let vertices = vec![
    ///     Vector::ZERO,
    ///     Vector::new(1.0, 1.0, 0.0),
    ///     Vector::new(2.0, 0.0, 0.0),
    ///     Vector::new(3.0, 1.0, 0.0),
    /// ];
    /// let polyline = Polyline::new(vertices, None);
    /// assert_eq!(polyline.num_segments(), 3);
    /// # }
    /// ```
    ///
    /// # Custom Connectivity Example
    ///
    /// ```
    /// # #[cfg(all(feature = "dim2", feature = "f32"))] {
    /// use parry2d::shape::Polyline;
    /// use parry2d::math::Vector;
    ///
    /// // Create a square with custom indices
    /// let vertices = vec![
    ///     Vector::ZERO,
    ///     Vector::new(1.0, 0.0),
    ///     Vector::new(1.0, 1.0),
    ///     Vector::new(0.0, 1.0),
    /// ];
    ///
    /// // Define edges to form a closed square
    /// let indices = vec![
    ///     [0, 1], [1, 2], [2, 3], [3, 0]
    /// ];
    ///
    /// let square = Polyline::new(vertices, Some(indices));
    /// assert_eq!(square.num_segments(), 4);
    ///
    /// // Each segment connects the correct vertices
    /// let first_segment = square.segment(0);
    /// assert_eq!(first_segment.a, Vector::ZERO);
    /// assert_eq!(first_segment.b, Vector::new(1.0, 0.0));
    /// # }
    /// ```
    pub fn new(vertices: Vec<Vector>, indices: Option<Vec<[u32; 2]>>) -> Self {
        // Index from 1 so empty input produces no segments instead of underflowing on `len - 1`.
        let indices =
            indices.unwrap_or_else(|| (1..vertices.len() as u32).map(|i| [i - 1, i]).collect());
        let leaves = indices.iter().enumerate().map(|(i, idx)| {
            let aabb =
                Segment::new(vertices[idx[0] as usize], vertices[idx[1] as usize]).local_aabb();
            (i, aabb)
        });

        // NOTE: we apply no dilation factor because we won't
        // update this tree dynamically.
        let bvh = Bvh::from_iter(BvhBuildStrategy::Binned, leaves);

        Self {
            bvh,
            vertices,
            indices,
            #[cfg(feature = "dim2")]
            pseudo_normals: None,
            #[cfg(feature = "dim2")]
            flags: PolylineFlags::empty(),
        }
    }

    /// Creates a new polyline with the given [`PolylineFlags`] controlling its optional associated
    /// data, e.g. orientation via [`PolylineFlags::ORIENTED`].
    ///
    /// # Example
    ///
    /// ```
    /// # #[cfg(all(feature = "dim2", feature = "f32"))] {
    /// use parry2d::shape::{Polyline, PolylineFlags};
    /// use parry2d::math::Vector;
    ///
    /// // A unit square wound counter-clockwise, so the solid is inside and outward points away.
    /// let vertices = vec![
    ///     Vector::new(-1.0, -1.0),
    ///     Vector::new(1.0, -1.0),
    ///     Vector::new(1.0, 1.0),
    ///     Vector::new(-1.0, 1.0),
    /// ];
    /// let indices = vec![[0, 1], [1, 2], [2, 3], [3, 0]];
    /// let polyline = Polyline::with_flags(vertices, Some(indices), PolylineFlags::ORIENTED);
    ///
    /// // The bottom edge's outward normal points down, away from the interior.
    /// let bottom = polyline.segment_normal_constraints(0).unwrap();
    /// assert!(bottom.face.abs_diff_eq(Vector::new(0.0, -1.0), 1.0e-5));
    /// # }
    /// ```
    #[cfg(feature = "dim2")]
    pub fn with_flags(
        vertices: Vec<Vector>,
        indices: Option<Vec<[u32; 2]>>,
        flags: PolylineFlags,
    ) -> Self {
        let mut result = Self::new(vertices, indices);
        result.set_flags(flags);
        result
    }

    /// Sets the [`PolylineFlags`], computing or discarding the polyline's optional associated data.
    #[cfg(feature = "dim2")]
    pub fn set_flags(&mut self, flags: PolylineFlags) {
        self.flags = flags;

        if flags.contains(PolylineFlags::ORIENTED) {
            self.compute_pseudo_normals();
        } else {
            self.pseudo_normals = None;
        }
    }

    /// The [`PolylineFlags`] controlling this polyline's optional associated data.
    #[cfg(feature = "dim2")]
    pub fn flags(&self) -> PolylineFlags {
        self.flags
    }

    /// Computes the outward pseudo-normal at every vertex (the normalized sum of its incident
    /// segments' outward normals) for the one-sided behavior of [`PolylineFlags::ORIENTED`].
    #[cfg(feature = "dim2")]
    fn compute_pseudo_normals(&mut self) {
        let mut vertex_normals = Vec::new();
        vertex_normals.resize(self.vertices.len(), Vector::ZERO);

        // A 2D vertex has at most two incident segments, so the normalized sum is their exact
        // bisector -- no angle weighting (unlike the 3D `TrianglePseudoNormals`).
        for idx in &self.indices {
            let a = idx[0] as usize;
            let b = idx[1] as usize;
            let normal = crate::utils::ccw_face_normal([self.vertices[a], self.vertices[b]])
                .unwrap_or(Vector::ZERO);
            vertex_normals[a] += normal;
            vertex_normals[b] += normal;
        }

        for normal in &mut vertex_normals {
            *normal = normal.normalize_or_zero();
        }

        self.pseudo_normals = Some(vertex_normals);
    }

    /// Returns the [`SegmentPseudoNormals`] for the segment with index `i`, or `None` unless this
    /// polyline was built with [`PolylineFlags::ORIENTED`].
    ///
    /// # Example
    ///
    /// ```
    /// # #[cfg(all(feature = "dim2", feature = "f32"))] {
    /// use parry2d::shape::{Polyline, PolylineFlags};
    /// use parry2d::math::Vector;
    ///
    /// let vertices = vec![Vector::new(0.0, 0.0), Vector::new(2.0, 0.0), Vector::new(2.0, 2.0)];
    /// let mut polyline = Polyline::new(vertices, Some(vec![[0, 1], [1, 2]]));
    ///
    /// // A double-sided polyline has no constraints...
    /// assert!(polyline.segment_normal_constraints(0).is_none());
    ///
    /// // ...until it is oriented.
    /// polyline.set_flags(PolylineFlags::ORIENTED);
    /// assert!(polyline.segment_normal_constraints(0).is_some());
    /// # }
    /// ```
    #[cfg(feature = "dim2")]
    pub fn segment_normal_constraints(&self, i: u32) -> Option<SegmentPseudoNormals> {
        let pseudo_normals = self.pseudo_normals.as_ref()?;
        let idx = self.indices[i as usize];
        let a = idx[0] as usize;
        let b = idx[1] as usize;
        let face = crate::utils::ccw_face_normal([self.vertices[a], self.vertices[b]])?;
        Some(SegmentPseudoNormals {
            face,
            edges: [pseudo_normals[a], pseudo_normals[b]],
        })
    }

    /// Pseudo-normals are a 2D-only feature; in 3D there are no segment normal constraints. 3D stub
    /// so the composite-shape impls compile; mirrors `TriMesh::triangle_normal_constraints`'s dim2
    /// stub.
    #[cfg(feature = "dim3")]
    #[doc(hidden)]
    pub fn segment_normal_constraints(&self, _i: u32) -> Option<SegmentPseudoNormals> {
        None
    }

    /// Computes the axis-aligned bounding box of this polyline in world space.
    ///
    /// The AABB is the smallest box aligned with the world axes that fully contains
    /// the polyline after applying the given position/rotation transformation.
    ///
    /// # Arguments
    ///
    /// * `pos` - The position and orientation (isometry) of the polyline in world space
    ///
    /// # Returns
    ///
    /// An `Aabb` that bounds the transformed polyline
    ///
    /// # Example
    ///
    /// ```
    /// # #[cfg(all(feature = "dim3", feature = "f32"))] {
    /// use parry3d::shape::Polyline;
    /// use parry3d::math::{Vector, Pose};
    ///
    /// // Create a polyline along the X axis
    /// let vertices = vec![
    ///     Vector::ZERO,
    ///     Vector::new(2.0, 0.0, 0.0),
    /// ];
    /// let polyline = Polyline::new(vertices, None);
    ///
    /// // Compute AABB at the origin
    /// let identity = Pose::identity();
    /// let aabb = polyline.aabb(&identity);
    /// assert_eq!(aabb.mins.x, 0.0);
    /// assert_eq!(aabb.maxs.x, 2.0);
    ///
    /// // Compute AABB after translating by (10, 5, 0)
    /// let translated = Pose::translation(10.0, 5.0, 0.0);
    /// let aabb_translated = polyline.aabb(&translated);
    /// assert_eq!(aabb_translated.mins.x, 10.0);
    /// assert_eq!(aabb_translated.maxs.x, 12.0);
    /// # }
    /// ```
    pub fn aabb(&self, pos: &Pose) -> Aabb {
        self.bvh.root_aabb().transform_by(pos)
    }

    /// Gets the local axis-aligned bounding box of this polyline.
    ///
    /// This returns the AABB in the polyline's local coordinate system (before any
    /// transformation is applied). It's more efficient than `aabb()` when you don't
    /// need to transform the polyline.
    ///
    /// # Returns
    ///
    /// An `Aabb` that bounds the polyline in local space
    ///
    /// # Example
    ///
    /// ```
    /// # #[cfg(all(feature = "dim2", feature = "f32"))] {
    /// use parry2d::shape::Polyline;
    /// use parry2d::math::Vector;
    ///
    /// // Create a rectangular polyline
    /// let vertices = vec![
    ///     Vector::new(-1.0, -2.0),
    ///     Vector::new(3.0, -2.0),
    ///     Vector::new(3.0, 4.0),
    ///     Vector::new(-1.0, 4.0),
    /// ];
    /// let polyline = Polyline::new(vertices, None);
    ///
    /// // Get the local AABB
    /// let aabb = polyline.local_aabb();
    ///
    /// // The AABB should contain all vertices
    /// assert_eq!(aabb.mins.x, -1.0);
    /// assert_eq!(aabb.mins.y, -2.0);
    /// assert_eq!(aabb.maxs.x, 3.0);
    /// assert_eq!(aabb.maxs.y, 4.0);
    /// # }
    /// ```
    pub fn local_aabb(&self) -> Aabb {
        self.bvh.root_aabb()
    }

    /// The BVH acceleration structure for this polyline.
    pub fn bvh(&self) -> &Bvh {
        &self.bvh
    }

    /// Returns the number of segments (edges) in this polyline.
    ///
    /// Each segment connects two vertices. For a polyline with `n` vertices and
    /// sequential connectivity, there are `n-1` segments. For custom connectivity,
    /// the number of segments equals the number of index pairs.
    ///
    /// # Returns
    ///
    /// The total number of line segments
    ///
    /// # Example
    ///
    /// ```
    /// # #[cfg(all(feature = "dim3", feature = "f32"))] {
    /// use parry3d::shape::Polyline;
    /// use parry3d::math::Vector;
    ///
    /// // Sequential polyline: 5 vertices -> 4 segments
    /// let vertices = vec![
    ///     Vector::ZERO,
    ///     Vector::new(1.0, 0.0, 0.0),
    ///     Vector::new(2.0, 0.0, 0.0),
    ///     Vector::new(3.0, 0.0, 0.0),
    ///     Vector::new(4.0, 0.0, 0.0),
    /// ];
    /// let polyline = Polyline::new(vertices.clone(), None);
    /// assert_eq!(polyline.num_segments(), 4);
    ///
    /// // Custom connectivity: can have different number of segments
    /// let indices = vec![[0, 4], [1, 3]]; // Only 2 segments
    /// let custom = Polyline::new(vertices, Some(indices));
    /// assert_eq!(custom.num_segments(), 2);
    /// # }
    /// ```
    pub fn num_segments(&self) -> usize {
        self.indices.len()
    }

    /// Returns an iterator over all segments in this polyline.
    ///
    /// Each segment is returned as a [`Segment`] object with two endpoints.
    /// The iterator yields exactly `num_segments()` items.
    ///
    /// # Returns
    ///
    /// An exact-size iterator that yields `Segment` instances
    ///
    /// # Example
    ///
    /// ```
    /// # #[cfg(all(feature = "dim2", feature = "f32"))] {
    /// use parry2d::shape::Polyline;
    /// use parry2d::math::Vector;
    ///
    /// // Create a triangle
    /// let vertices = vec![
    ///     Vector::ZERO,
    ///     Vector::new(1.0, 0.0),
    ///     Vector::new(0.5, 1.0),
    /// ];
    /// let polyline = Polyline::new(vertices, None);
    ///
    /// // Iterate over all segments
    /// let mut total_length = 0.0;
    /// for segment in polyline.segments() {
    ///     total_length += segment.length();
    /// }
    ///
    /// // Calculate expected perimeter (not closed, so 2 sides only)
    /// assert!(total_length > 2.0);
    ///
    /// // Collect all segments into a vector
    /// let segments: Vec<_> = polyline.segments().collect();
    /// assert_eq!(segments.len(), 2);
    /// # }
    /// ```
    pub fn segments(&self) -> impl ExactSizeIterator<Item = Segment> + '_ {
        self.indices.iter().map(move |ids| {
            Segment::new(
                self.vertices[ids[0] as usize],
                self.vertices[ids[1] as usize],
            )
        })
    }

    /// Returns the segment at the given index.
    ///
    /// This retrieves a specific segment by its index. Indices range from `0` to
    /// `num_segments() - 1`. If you need to access multiple segments, consider
    /// using the `segments()` iterator instead.
    ///
    /// # Arguments
    ///
    /// * `i` - The index of the segment to retrieve (0-based)
    ///
    /// # Returns
    ///
    /// A `Segment` representing the edge at index `i`
    ///
    /// # Panics
    ///
    /// Panics if `i >= num_segments()`
    ///
    /// # Example
    ///
    /// ```
    /// # #[cfg(all(feature = "dim3", feature = "f32"))] {
    /// use parry3d::shape::Polyline;
    /// use parry3d::math::Vector;
    ///
    /// let vertices = vec![
    ///     Vector::ZERO,
    ///     Vector::new(1.0, 0.0, 0.0),
    ///     Vector::new(2.0, 1.0, 0.0),
    /// ];
    /// let polyline = Polyline::new(vertices, None);
    ///
    /// // Get the first segment (connects vertex 0 to vertex 1)
    /// let seg0 = polyline.segment(0);
    /// assert_eq!(seg0.a, Vector::ZERO);
    /// assert_eq!(seg0.b, Vector::new(1.0, 0.0, 0.0));
    /// assert_eq!(seg0.length(), 1.0);
    ///
    /// // Get the second segment (connects vertex 1 to vertex 2)
    /// let seg1 = polyline.segment(1);
    /// assert_eq!(seg1.a, Vector::new(1.0, 0.0, 0.0));
    /// assert_eq!(seg1.b, Vector::new(2.0, 1.0, 0.0));
    /// # }
    /// ```
    pub fn segment(&self, i: u32) -> Segment {
        let idx = self.indices[i as usize];
        Segment::new(
            self.vertices[idx[0] as usize],
            self.vertices[idx[1] as usize],
        )
    }

    /// Transforms  the feature-id of a segment to the feature-id of this polyline.
    pub fn segment_feature_to_polyline_feature(
        &self,
        segment: u32,
        _feature: FeatureId,
    ) -> FeatureId {
        // TODO: return a vertex feature when it makes sense.
        #[cfg(feature = "dim2")]
        return FeatureId::Face(segment);
        #[cfg(feature = "dim3")]
        return FeatureId::Edge(segment);
    }

    /// Returns a slice containing all vertices of this polyline.
    ///
    /// Vertices are the points that define the polyline. Segments connect
    /// pairs of these vertices according to the index buffer.
    ///
    /// # Returns
    ///
    /// A slice of all vertex points
    ///
    /// # Example
    ///
    /// ```
    /// # #[cfg(all(feature = "dim2", feature = "f32"))] {
    /// use parry2d::shape::Polyline;
    /// use parry2d::math::Vector;
    ///
    /// let vertices = vec![
    ///     Vector::ZERO,
    ///     Vector::new(1.0, 0.0),
    ///     Vector::new(1.0, 1.0),
    /// ];
    /// let polyline = Polyline::new(vertices.clone(), None);
    ///
    /// // Access all vertices
    /// let verts = polyline.vertices();
    /// assert_eq!(verts.len(), 3);
    /// assert_eq!(verts[0], Vector::ZERO);
    /// assert_eq!(verts[1], Vector::new(1.0, 0.0));
    /// assert_eq!(verts[2], Vector::new(1.0, 1.0));
    ///
    /// // You can iterate over vertices
    /// for (i, vertex) in polyline.vertices().iter().enumerate() {
    ///     println!("Vertex {}: {:?}", i, vertex);
    /// }
    /// # }
    /// ```
    pub fn vertices(&self) -> &[Vector] {
        &self.vertices[..]
    }

    /// Returns a slice containing all segment indices.
    ///
    /// Each index is a pair `[u32; 2]` representing a segment connecting two vertices.
    /// The first element is the index of the segment's start vertex, and the second
    /// is the index of the end vertex.
    ///
    /// # Returns
    ///
    /// A slice of all segment index pairs
    ///
    /// # Example
    ///
    /// ```
    /// # #[cfg(all(feature = "dim3", feature = "f32"))] {
    /// use parry3d::shape::Polyline;
    /// use parry3d::math::Vector;
    ///
    /// let vertices = vec![
    ///     Vector::ZERO,
    ///     Vector::new(1.0, 0.0, 0.0),
    ///     Vector::new(2.0, 0.0, 0.0),
    /// ];
    ///
    /// // With automatic indices
    /// let polyline = Polyline::new(vertices.clone(), None);
    /// let indices = polyline.indices();
    /// assert_eq!(indices.len(), 2);
    /// assert_eq!(indices[0], [0, 1]); // First segment: vertex 0 -> 1
    /// assert_eq!(indices[1], [1, 2]); // Second segment: vertex 1 -> 2
    ///
    /// // With custom indices
    /// let custom_indices = vec![[0, 2], [1, 0]];
    /// let custom = Polyline::new(vertices, Some(custom_indices));
    /// assert_eq!(custom.indices()[0], [0, 2]);
    /// assert_eq!(custom.indices()[1], [1, 0]);
    /// # }
    /// ```
    pub fn indices(&self) -> &[[u32; 2]] {
        &self.indices
    }

    /// A flat view of the index buffer of this mesh.
    pub fn flat_indices(&self) -> &[u32] {
        unsafe {
            let len = self.indices.len() * 2;
            let data = self.indices.as_ptr() as *const u32;
            core::slice::from_raw_parts(data, len)
        }
    }

    /// Computes a scaled version of this polyline.
    ///
    /// This consumes the polyline and returns a new one with all vertices scaled
    /// component-wise by the given scale vector. The connectivity (indices) remains
    /// unchanged, but the BVH is rebuilt to reflect the new geometry.
    ///
    /// # Arguments
    ///
    /// * `scale` - The scaling factors for each axis
    ///
    /// # Returns
    ///
    /// A new polyline with scaled vertices
    ///
    /// # Example
    ///
    /// ```
    /// # #[cfg(all(feature = "dim2", feature = "f32"))] {
    /// use parry2d::shape::Polyline;
    /// use parry2d::math::Vector;
    ///
    /// let vertices = vec![
    ///     Vector::new(1.0, 2.0),
    ///     Vector::new(3.0, 4.0),
    /// ];
    /// let polyline = Polyline::new(vertices, None);
    ///
    /// // Scale by 2x in X and 3x in Y
    /// let scaled = polyline.scaled(Vector::new(2.0, 3.0));
    ///
    /// // Check scaled vertices
    /// assert_eq!(scaled.vertices()[0], Vector::new(2.0, 6.0));
    /// assert_eq!(scaled.vertices()[1], Vector::new(6.0, 12.0));
    /// # }
    /// ```
    ///
    /// # Note
    ///
    /// This method consumes `self`. If you need to keep the original polyline,
    /// clone it first:
    ///
    /// ```
    /// # #[cfg(all(feature = "dim3", feature = "f32"))] {
    /// use parry3d::shape::Polyline;
    /// use parry3d::math::Vector;
    ///
    /// let vertices = vec![Vector::ZERO, Vector::new(1.0, 0.0, 0.0)];
    /// let original = Polyline::new(vertices, None);
    ///
    /// // Clone before scaling if you need to keep the original
    /// let scaled = original.clone().scaled(Vector::new(2.0, 2.0, 2.0));
    ///
    /// // Both polylines still exist
    /// assert_eq!(original.vertices()[1].x, 1.0);
    /// assert_eq!(scaled.vertices()[1].x, 2.0);
    /// # }
    /// ```
    pub fn scaled(mut self, scale: Vector) -> Self {
        self.vertices.iter_mut().for_each(|pt| *pt *= scale);
        let mut bvh = self.bvh.clone();
        bvh.scale(scale);

        #[cfg(feature = "dim2")]
        {
            let mut result = Self {
                bvh,
                vertices: self.vertices,
                indices: self.indices,
                pseudo_normals: None,
                flags: PolylineFlags::empty(),
            };
            // Recompute the pseudo-normals from the scaled geometry.
            result.set_flags(self.flags);
            result
        }
        #[cfg(feature = "dim3")]
        {
            Self {
                bvh,
                vertices: self.vertices,
                indices: self.indices,
            }
        }
    }

    /// Reverses the orientation of this polyline.
    ///
    /// This operation:
    /// 1. Swaps the start and end vertex of each segment (reversing edge direction)
    /// 2. Reverses the order of segments in the index buffer
    /// 3. Rebuilds the BVH to maintain correct acceleration structure
    ///
    /// After reversing, traversing the polyline segments in order will visit
    /// the same geometry but in the opposite direction.
    ///
    /// # Example
    ///
    /// ```
    /// # #[cfg(all(feature = "dim2", feature = "f32"))] {
    /// use parry2d::shape::Polyline;
    /// use parry2d::math::Vector;
    ///
    /// let vertices = vec![
    ///     Vector::ZERO,
    ///     Vector::new(1.0, 0.0),
    ///     Vector::new(2.0, 0.0),
    /// ];
    /// let mut polyline = Polyline::new(vertices, None);
    ///
    /// // Original: segment 0 goes from vertex 0 to 1, segment 1 from 1 to 2
    /// assert_eq!(polyline.indices()[0], [0, 1]);
    /// assert_eq!(polyline.indices()[1], [1, 2]);
    ///
    /// // Reverse the polyline
    /// polyline.reverse();
    ///
    /// // After reversing: order is flipped and directions are swapped
    /// // The last segment becomes first, with swapped endpoints
    /// assert_eq!(polyline.indices()[0], [2, 1]);
    /// assert_eq!(polyline.indices()[1], [1, 0]);
    /// # }
    /// ```
    ///
    /// # Use Cases
    ///
    /// This is useful for:
    /// - Correcting winding order for 2D shapes
    /// - Reversing path direction for navigation
    /// - Ensuring consistent edge orientation in connected components
    pub fn reverse(&mut self) {
        for idx in &mut self.indices {
            idx.swap(0, 1);
        }

        self.indices.reverse();

        // Rebuild the bvh since the segment indices no longer map to the correct element.
        // TODO PERF: should the Bvh have a function for efficient leaf index remapping?
        //            Probably not worth it unless this function starts showing up as a
        //            bottleneck for someone.
        let leaves = self.segments().map(|seg| seg.local_aabb()).enumerate();
        let bvh = Bvh::from_iter(BvhBuildStrategy::Binned, leaves);
        self.bvh = bvh;

        // Reversing flips the winding, so recompute the outward side.
        #[cfg(feature = "dim2")]
        if self.flags.contains(PolylineFlags::ORIENTED) {
            self.compute_pseudo_normals();
        }
    }

    /// Extracts the connected components of this polyline, consuming `self`.
    ///
    /// This method is currently quite restrictive on the kind of allowed input. The polyline
    /// represented by `self` must already have an index buffer sorted such that:
    /// - Each connected component appears in the index buffer one after the other, i.e., a
    ///   connected component of this polyline must be a contiguous range of this polyline’s
    ///   index buffer.
    /// - Each connected component is closed, i.e., each range of this polyline index buffer
    ///   `self.indices[i_start..=i_end]` forming a complete connected component, we must have
    ///   `self.indices[i_start][0] == self.indices[i_end][1]`.
    /// - The indices for each component must already be in order, i.e., if the segments
    ///   `self.indices[i]` and `self.indices[i + 1]` are part of the same connected component then
    ///   we must have `self.indices[i][1] == self.indices[i + 1][0]`.
    ///
    /// # Output
    /// Returns the set of polylines. If the inputs fulfill the constraints mentioned above, each
    /// polyline will be a closed loop with consistent edge orientations, i.e., for all indices `i`,
    /// we have `polyline.indices[i][1] == polyline.indices[i + 1][0]`.
    ///
    /// The orientation of each closed loop (clockwise or counterclockwise) are identical to their
    /// original orientation in `self`.
    pub fn extract_connected_components(&self) -> Vec<Polyline> {
        let vertices = self.vertices();
        let indices = self.indices();

        if indices.is_empty() {
            // Polyline is empty, return empty Vec
            Vec::new()
        } else {
            let mut components = Vec::new();

            let mut start_i = 0; // Start position of component
            let mut start_node = indices[0][0]; // Start vertex index of component

            let mut component_vertices = Vec::new();
            let mut component_indices: Vec<[u32; 2]> = Vec::new();

            // Iterate over indices, building polylines as we go
            for (i, idx) in indices.iter().enumerate() {
                component_vertices.push(vertices[idx[0] as usize]);

                if idx[1] != start_node {
                    // Keep scanning and adding data
                    component_indices.push([(i - start_i) as u32, (i - start_i + 1) as u32]);
                } else {
                    // Start node reached: build polyline and start next component
                    component_indices.push([(i - start_i) as u32, 0]);
                    components.push(Polyline::new(
                        core::mem::take(&mut component_vertices),
                        Some(core::mem::take(&mut component_indices)),
                    ));

                    if i + 1 < indices.len() {
                        // More components to find
                        start_node = indices[i + 1][0];
                        start_i = i + 1;
                    }
                }
            }

            components
        }
    }

    /// Perform a point projection assuming a solid interior based on a counter-clock-wise orientation.
    ///
    /// This is similar to `self.project_local_point_and_get_location` except that the resulting
    /// `PointProjection::is_inside` will be set to true if the point is inside of the area delimited
    /// by this polyline, assuming that:
    /// - This polyline isn’t self-crossing.
    /// - This polyline is closed with `self.indices[i][1] == self.indices[(i + 1) % num_indices][0]` where
    ///   `num_indices == self.indices.len()`.
    /// - This polyline is oriented counter-clockwise.
    /// - In 3D, the polyline is assumed to be fully coplanar, on a plane with normal given by
    ///   `axis`.
    ///
    /// These properties are not checked.
    pub fn project_local_point_assuming_solid_interior_ccw(
        &self,
        point: Vector,
        #[cfg(feature = "dim3")] axis: u8,
    ) -> (PointProjection, (u32, SegmentPointLocation)) {
        let mut proj = self.project_local_point_and_get_location(point, false);
        let segment1 = self.segment((proj.1).0);

        #[cfg(feature = "dim2")]
        let normal1 = segment1.normal();
        #[cfg(feature = "dim3")]
        let normal1 = segment1.planar_normal(axis);

        if let Some(normal1) = normal1 {
            proj.0.is_inside = match proj.1 .1 {
                SegmentPointLocation::OnVertex(i) => {
                    let dir2 = if i == 0 {
                        let adj_seg = if proj.1 .0 == 0 {
                            self.indices().len() as u32 - 1
                        } else {
                            proj.1 .0 - 1
                        };

                        assert_eq!(segment1.a, self.segment(adj_seg).b);
                        -self.segment(adj_seg).scaled_direction()
                    } else {
                        assert_eq!(i, 1);
                        let adj_seg = (proj.1 .0 + 1) % self.indices().len() as u32;
                        assert_eq!(segment1.b, self.segment(adj_seg).a);

                        self.segment(adj_seg).scaled_direction()
                    };

                    let dot = normal1.dot(dir2);
                    // TODO: is this threshold too big? This corresponds to an angle equal to
                    //       abs(acos(1.0e-3)) = (90 - 0.057) degrees.
                    //       We did encounter some cases where this was needed, but perhaps the
                    //       actual problem was an issue with the SegmentPointLocation (which should
                    //       perhaps have been Edge instead of Vertex)?
                    let threshold = 1.0e-3 * dir2.length();
                    if dot.abs() > threshold {
                        // If the vertex is a reentrant vertex, then the point is
                        // inside. Otherwise, it is outside.
                        dot >= 0.0
                    } else {
                        // If the two edges are collinear, we can’t classify the vertex.
                        // So check against the edge’s normal instead.
                        (point - proj.0.point).dot(normal1) <= 0.0
                    }
                }
                SegmentPointLocation::OnEdge(_) => (point - proj.0.point).dot(normal1) <= 0.0,
            };
        }

        proj
    }
}

impl CompositeShape for Polyline {
    fn map_part_at(
        &self,
        i: u32,
        f: &mut dyn FnMut(Option<&Pose>, &dyn Shape, Option<&dyn NormalConstraints>),
    ) {
        let seg = self.segment(i);
        let normals = self.segment_normal_constraints(i);
        f(
            None,
            &seg,
            normals.as_ref().map(|n| n as &dyn NormalConstraints),
        )
    }

    fn bvh(&self) -> &Bvh {
        &self.bvh
    }
}

impl TypedCompositeShape for Polyline {
    type PartShape = Segment;
    type PartNormalConstraints = SegmentPseudoNormals;

    #[inline(always)]
    fn map_typed_part_at<T>(
        &self,
        i: u32,
        mut f: impl FnMut(Option<&Pose>, &Self::PartShape, Option<&Self::PartNormalConstraints>) -> T,
    ) -> Option<T> {
        let seg = self.segment(i);
        let normals = self.segment_normal_constraints(i);
        Some(f(None, &seg, normals.as_ref()))
    }

    #[inline(always)]
    fn map_untyped_part_at<T>(
        &self,
        i: u32,
        mut f: impl FnMut(Option<&Pose>, &dyn Shape, Option<&dyn NormalConstraints>) -> T,
    ) -> Option<T> {
        let seg = self.segment(i);
        let normals = self.segment_normal_constraints(i);
        Some(f(
            None,
            &seg,
            normals.as_ref().map(|n| n as &dyn NormalConstraints),
        ))
    }
}

#[cfg(test)]
#[cfg(all(feature = "dim2", feature = "alloc"))]
mod pseudo_normal_tests {
    use crate::math::Vector;
    use crate::shape::{Polyline, PolylineFlags};

    fn ccw_square() -> Polyline {
        // CCW unit square: solid inside, outward away from the center.
        let vertices = vec![
            Vector::new(-1.0, -1.0),
            Vector::new(1.0, -1.0),
            Vector::new(1.0, 1.0),
            Vector::new(-1.0, 1.0),
        ];
        Polyline::new(vertices, Some(vec![[0, 1], [1, 2], [2, 3], [3, 0]]))
    }

    #[test]
    fn not_oriented_by_default() {
        assert!(ccw_square().segment_normal_constraints(0).is_none());
    }

    #[test]
    fn face_and_edge_normals_are_unit_and_outward() {
        let mut polyline = ccw_square();
        polyline.set_flags(PolylineFlags::ORIENTED);

        for i in 0..polyline.num_segments() as u32 {
            let segment = polyline.segment(i);
            let constraints = polyline.segment_normal_constraints(i).unwrap();

            assert!((constraints.face.length() - 1.0).abs() < 1.0e-5);
            for edge in constraints.edges {
                assert!((edge.length() - 1.0).abs() < 1.0e-5);
            }

            // The square is centered on the origin, so a midpoint doubles as its outward direction.
            let midpoint = (segment.a + segment.b) * 0.5;
            assert!(constraints.face.dot(midpoint) > 0.0);
        }
    }

    #[test]
    fn corner_pseudo_normal_bisects_its_two_faces() {
        let mut polyline = ccw_square();
        polyline.set_flags(PolylineFlags::ORIENTED);

        // Vertex 1's pseudo-normal bisects the bottom edge (-Y) and right edge (+X): (1, -1) normalized.
        let bottom = polyline.segment_normal_constraints(0).unwrap();
        assert!(bottom.edges[1].abs_diff_eq(Vector::new(1.0, -1.0).normalize(), 1.0e-5));
    }

    #[test]
    fn degenerate_input_yields_no_segments() {
        // Auto-generated indices must not underflow `len - 1` on 0- or 1-vertex input.
        assert_eq!(Polyline::new(vec![], None).num_segments(), 0);
        assert_eq!(Polyline::new(vec![Vector::ZERO], None).num_segments(), 0);
    }

    #[test]
    fn oriented_point_query_treats_interior_as_inside() {
        use crate::query::{PointQuery, PointQueryWithLocation};

        let mut polyline = ccw_square();
        polyline.set_flags(PolylineFlags::ORIENTED);
        let inside = Vector::new(0.25, -0.5);

        // Solid: an inside point is its own projection.
        let solid = polyline.project_local_point(inside, true);
        assert!(solid.is_inside);
        assert!(solid.point.abs_diff_eq(inside, 1.0e-5));

        // Non-solid: still inside, but projected to the boundary.
        let hollow = polyline.project_local_point(inside, false);
        assert!(hollow.is_inside);
        assert!(!hollow.point.abs_diff_eq(inside, 1.0e-5));

        assert!(polyline.contains_local_point(inside));
        assert!(
            polyline
                .project_local_point_and_get_location(inside, false)
                .0
                .is_inside
        );
        assert!(
            polyline
                .project_local_point_and_get_feature(inside)
                .0
                .is_inside
        );
    }

    #[test]
    fn oriented_point_query_leaves_exterior_outside() {
        use crate::query::PointQuery;

        let mut polyline = ccw_square();
        polyline.set_flags(PolylineFlags::ORIENTED);
        let outside = Vector::new(2.0, 0.0);

        let proj = polyline.project_local_point(outside, true);
        assert!(!proj.is_inside);
        assert!(proj.point.abs_diff_eq(Vector::new(1.0, 0.0), 1.0e-5));
        assert!(!polyline.contains_local_point(outside));
    }

    #[test]
    fn unoriented_polyline_has_no_interior() {
        use crate::query::PointQuery;

        // Unoriented: a hollow wireframe with no interior.
        let polyline = ccw_square();
        let center = Vector::ZERO;

        assert!(!polyline.project_local_point(center, true).is_inside);
        assert!(!polyline.contains_local_point(center));
    }
}
