//!
//! Shape composed from the union of primitives.
//!

use crate::bounding_volume::{Aabb, BoundingSphere, BoundingVolume};
use crate::math::Pose;
use crate::math::{Real, Vector};
use crate::partitioning::{Bvh, BvhBuildStrategy};
use crate::query::details::NormalConstraints;
#[cfg(feature = "dim2")]
use crate::shape::CompoundEdgeCone;
#[cfg(feature = "dim3")]
use crate::shape::CompoundFaceCone;
use crate::shape::CompoundPseudoNormals;
use crate::shape::{CompositeShape, Shape, SharedShape, TypedCompositeShape};
#[cfg(feature = "dim2")]
use crate::shape::{ConvexPolygon, TriMesh, Triangle};
#[cfg(feature = "dim2")]
use crate::transformation::hertel_mehlhorn;
use crate::utils::hashmap::HashMap;
use alloc::vec::Vec;

#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[repr(C)]
#[derive(Clone, Copy, Debug, Default, Eq, Hash, Ord, PartialEq, PartialOrd)]
/// Controls how a [`Compound`] is loaded.
pub struct CompoundFlags(u8);

bitflags::bitflags! {
    impl CompoundFlags: u8 {
        /// If set, the edges (2D) or faces (3D) where two parts meet are treated as interior to
        /// the union, and contact normals are clamped to the surviving outline. This removes the
        /// ledge a body catches on when it slides across the cut between two parts of a convex
        /// decomposition.
        ///
        /// This is the [`Compound`] counterpart of `TriMeshFlags::FIX_INTERNAL_EDGES`, and like it
        /// costs a one-off pass over the parts' edges when set.
        const FIX_INTERNAL_EDGES = 1;
    }
}

/// A compound shape with an aabb bounding volume.
///
/// A compound shape is a shape composed of the union of several simpler shape. This is
/// the main way of creating a concave shape from convex parts. Each parts can have its own
/// delta transformation to shift or rotate it with regard to the other shapes.
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[derive(Clone, Debug)]
pub struct Compound {
    shapes: Vec<(Pose, SharedShape)>,
    bvh: Bvh,
    aabbs: Vec<Aabb>,
    aabb: Aabb,
    #[cfg_attr(feature = "serde", serde(default))]
    flags: CompoundFlags,
    /// `None` for a part that constrains nothing: one that is not polygonal, or any part at all
    /// while [`CompoundFlags::FIX_INTERNAL_EDGES`] is unset.
    #[cfg_attr(feature = "serde", serde(default))]
    pseudo_normals: Option<Vec<Option<CompoundPseudoNormals>>>,
}

impl Compound {
    /// Builds a new compound shape from a collection of sub-shapes.
    ///
    /// A compound shape combines multiple primitive shapes into a single composite shape,
    /// each with its own relative position and orientation. This is the primary way to
    /// create concave shapes from convex primitives. The compound internally builds a
    /// Bounding Volume Hierarchy (BVH) for efficient collision detection.
    ///
    /// # Arguments
    ///
    /// * `shapes` - A vector of (position, shape) pairs. Each pair defines:
    ///   - An [`Pose`] representing the sub-shape's position and orientation relative to the compound's origin
    ///   - A [`SharedShape`] containing the actual geometry
    ///
    /// # Panics
    ///
    /// - If the input vector is empty (a compound must contain at least one shape)
    /// - If any of the provided shapes are themselves composite shapes (nested composites are not allowed)
    ///
    /// # Performance
    ///
    /// The BVH is built using a binned construction strategy optimized for static scenes.
    /// For large compounds (100+ shapes), construction may take noticeable time but provides
    /// excellent query performance.
    ///
    /// # Example
    ///
    /// ```
    /// # #[cfg(all(feature = "dim3", feature = "f32"))] {
    /// # use parry3d::shape::{Compound, Ball, Cuboid, SharedShape};
    /// use parry3d::math::Pose;
    /// use parry3d::math::Vector;
    ///
    /// // Create a compound shape resembling a dumbbell
    /// let shapes = vec![
    ///     // Left sphere
    ///     (
    ///         Pose::translation(-2.0, 0.0, 0.0),
    ///         SharedShape::new(Ball::new(0.5))
    ///     ),
    ///     // Center bar
    ///     (
    ///         Pose::identity(),
    ///         SharedShape::new(Cuboid::new(Vector::new(2.0, 0.2, 0.2)))
    ///     ),
    ///     // Right sphere
    ///     (
    ///         Pose::translation(2.0, 0.0, 0.0),
    ///         SharedShape::new(Ball::new(0.5))
    ///     ),
    /// ];
    ///
    /// let compound = Compound::new(shapes);
    ///
    /// // The compound now contains all three shapes
    /// assert_eq!(compound.shapes().len(), 3);
    /// # }
    /// ```
    ///
    /// ```
    /// # #[cfg(all(feature = "dim2", feature = "f32"))] {
    /// # use parry2d::shape::{Compound, Ball, Cuboid, SharedShape};
    /// # use parry2d::math::Pose;
    /// # use parry2d::math::Vector;
    ///
    /// // Create an L-shaped compound
    /// let shapes = vec![
    ///     // Vertical rectangle
    ///     (
    ///         Pose::translation(0.0, 1.0),
    ///         SharedShape::new(Cuboid::new(Vector::new(0.5, 1.0)))
    ///     ),
    ///     // Horizontal rectangle
    ///     (
    ///         Pose::translation(1.0, 0.0),
    ///         SharedShape::new(Cuboid::new(Vector::new(1.0, 0.5)))
    ///     ),
    /// ];
    ///
    /// let l_shape = Compound::new(shapes);
    /// assert_eq!(l_shape.shapes().len(), 2);
    /// # }
    /// ```
    pub fn new(shapes: Vec<(Pose, SharedShape)>) -> Compound {
        assert!(
            !shapes.is_empty(),
            "A compound shape must contain at least one shape."
        );
        let mut aabbs = Vec::new();
        let mut leaves = Vec::new();
        let mut aabb = Aabb::new_invalid();

        for (i, (delta, shape)) in shapes.iter().enumerate() {
            let bv = shape.compute_aabb(delta);

            aabb.merge(&bv);
            aabbs.push(bv);
            leaves.push((i, bv));

            if shape.as_composite_shape().is_some() {
                panic!("Nested composite shapes are not allowed.");
            }
        }

        // NOTE: we apply no dilation factor because we won't
        // update this tree dynamically.
        let bvh = Bvh::from_iter(BvhBuildStrategy::Binned, leaves);

        Compound {
            shapes,
            bvh,
            aabbs,
            aabb,
            flags: CompoundFlags::empty(),
            pseudo_normals: None,
        }
    }

    /// Sets the [`CompoundFlags`], computing or discarding the compound's optional associated data.
    pub fn set_flags(&mut self, flags: CompoundFlags) {
        self.flags = flags;

        if flags.contains(CompoundFlags::FIX_INTERNAL_EDGES) {
            self.compute_pseudo_normals();
        } else {
            self.pseudo_normals = None;
        }
    }

    /// The [`CompoundFlags`] controlling this compound's optional associated data.
    pub fn flags(&self) -> CompoundFlags {
        self.flags
    }

    /// The parts' outlines in the compound's frame, wound counter-clockwise. `None` for a part
    /// with no straight sides, which therefore shares none.
    #[cfg(feature = "dim2")]
    fn part_outlines(&self) -> Vec<Option<Vec<Vector>>> {
        self.shapes
            .iter()
            .map(|(part_pos, part)| {
                let corners: Vec<Vector> = if let Some(polygon) = part.as_convex_polygon() {
                    polygon.points().to_vec()
                } else if let Some(cuboid) = part.as_cuboid() {
                    cuboid.to_polyline().to_vec()
                } else {
                    let triangle = part.as_triangle()?;
                    alloc::vec![triangle.a, triangle.b, triangle.c]
                };

                let mut outline: Vec<Vector> = corners
                    .into_iter()
                    .map(|corner| part_pos.transform_point(corner))
                    .collect();

                if signed_area(&outline) < 0.0 {
                    outline.reverse();
                }

                Some(outline)
            })
            .collect()
    }

    /// Builds the normal cones every part may report contacts in, from the outline of the union
    /// rather than from each part's own faces.
    ///
    /// An edge two parts share is a cut, so it is dropped, and each surviving edge opens only
    /// towards the neighbours that survived with it. That closes the corner where a surface runs
    /// into a cut -- the ledge a body would otherwise catch on.
    #[cfg(feature = "dim2")]
    fn compute_pseudo_normals(&mut self) {
        let outlines = self.part_outlines();

        // Parts meet on coordinates both were built from, so shared endpoints agree exactly and
        // the key can be the bits. Quantizing would merge edges that merely lie close together,
        // which is how a real surface goes missing. Adding zero collapses -0.0 onto 0.0.
        let edge_key = |a: Vector, b: Vector| {
            let (mut first, mut second) = (
                a.to_array()
                    .map(|coord| without_negative_zero(coord).to_bits()),
                b.to_array()
                    .map(|coord| without_negative_zero(coord).to_bits()),
            );
            if second < first {
                core::mem::swap(&mut first, &mut second);
            }
            (first, second)
        };

        let mut parts_sharing_edge = HashMap::default();
        for outline in outlines.iter().flatten() {
            for (a, b) in outline_edges(outline) {
                *parts_sharing_edge.entry(edge_key(a, b)).or_insert(0u32) += 1;
            }
        }

        let pseudo_normals = self
            .shapes
            .iter()
            .zip(outlines.iter())
            .map(|((part_pos, _), outline)| {
                let outline = outline.as_ref()?;

                let boundary_normals: Vec<Option<Vector>> = outline_edges(outline)
                    .map(|(a, b)| {
                        let shared = parts_sharing_edge
                            .get(&edge_key(a, b))
                            .is_some_and(|parts| *parts > 1);

                        (!shared)
                            .then(|| crate::utils::ccw_face_normal([a, b]))
                            .flatten()
                    })
                    .collect();

                let into_part = part_pos.rotation.inverse();
                let count = boundary_normals.len();
                let mut boundary_edges = Vec::new();

                for (i, face) in boundary_normals.iter().enumerate() {
                    let Some(face) = *face else {
                        continue;
                    };

                    // Halfway, so the two edges at a corner split its range and neither cone grows
                    // past a half turn, which the containment test could not represent.
                    let halfway_to = |neighbour: Option<Vector>| match neighbour {
                        Some(neighbour) => (face + neighbour).normalize_or(face),
                        None => face,
                    };

                    // Walking counter-clockwise turns the normals the same way, so the previous
                    // edge lies clockwise of this one. Cones are stated in the part's frame, which
                    // is where the narrow phase applies them.
                    boundary_edges.push(CompoundEdgeCone {
                        face: into_part * face,
                        clockwise_limit: into_part
                            * halfway_to(boundary_normals[(i + count - 1) % count]),
                        counter_clockwise_limit: into_part
                            * halfway_to(boundary_normals[(i + 1) % count]),
                    });
                }

                Some(CompoundPseudoNormals { boundary_edges })
            })
            .collect();

        self.pseudo_normals = Some(pseudo_normals);
    }

    /// The parts' faces in the compound's frame, each an outward normal with its vertex loop.
    /// `None` for a part that is not a polyhedron, which therefore shares no face.
    #[cfg(feature = "dim3")]
    #[allow(clippy::type_complexity)]
    fn part_faces(&self) -> Vec<Option<Vec<(Vector, Vec<Vector>)>>> {
        self.shapes
            .iter()
            .map(|(part_pos, part)| {
                let polyhedron = part.as_convex_polyhedron()?;
                let points = polyhedron.points();
                let corner_indices = polyhedron.vertices_adj_to_face();

                Some(
                    polyhedron
                        .faces()
                        .iter()
                        .map(|face| {
                            let first = face.first_vertex_or_edge as usize;
                            let corners = corner_indices
                                [first..first + face.num_vertices_or_edges as usize]
                                .iter()
                                .map(|index| part_pos.transform_point(points[*index as usize]))
                                .collect();
                            (part_pos.rotation * face.normal, corners)
                        })
                        .collect(),
                )
            })
            .collect()
    }

    /// Builds the normal cones every part may report contacts in, from the surface of the union
    /// rather than from each part's own faces.
    ///
    /// A face two parts share is a cut, so it is dropped, and each edge of a surviving face opens
    /// only towards the face that survived with it. That closes the corner where a surface runs
    /// into a cut -- the ledge a body would otherwise catch on.
    #[cfg(feature = "dim3")]
    fn compute_pseudo_normals(&mut self) {
        let part_faces = self.part_faces();

        // Parts of a decomposition meet on coordinates both were built from, so a shared face has
        // the same corners exactly and the key can be the bits, order ignored. Quantizing would
        // merge faces that only lie close together, which is how a real surface goes missing.
        let face_key = |corners: &[Vector]| {
            let mut key: Vec<_> = corners
                .iter()
                .map(|corner| {
                    corner
                        .to_array()
                        .map(|coord| without_negative_zero(coord).to_bits())
                })
                .collect();
            key.sort_unstable();
            key
        };
        let edge_key = |a: Vector, b: Vector| {
            let (mut first, mut second) = (
                a.to_array()
                    .map(|coord| without_negative_zero(coord).to_bits()),
                b.to_array()
                    .map(|coord| without_negative_zero(coord).to_bits()),
            );
            if second < first {
                core::mem::swap(&mut first, &mut second);
            }
            (first, second)
        };

        let mut parts_sharing_face = HashMap::default();
        for faces in part_faces.iter().flatten() {
            for (_, corners) in faces {
                *parts_sharing_face.entry(face_key(corners)).or_insert(0u32) += 1;
            }
        }

        let pseudo_normals = self
            .shapes
            .iter()
            .zip(part_faces.iter())
            .map(|((part_pos, _), faces)| {
                let faces = faces.as_ref()?;

                let is_boundary: Vec<bool> = faces
                    .iter()
                    .map(|(_, corners)| {
                        parts_sharing_face
                            .get(&face_key(corners))
                            .is_none_or(|parts| *parts <= 1)
                    })
                    .collect();

                // The boundary face on the other side of each edge, for the halfway pseudo-normal.
                let mut boundary_normal_across_edge = HashMap::default();
                for (face, boundary) in faces.iter().zip(is_boundary.iter()) {
                    if !boundary {
                        continue;
                    }
                    let (normal, corners) = face;
                    for i in 0..corners.len() {
                        let edge = edge_key(corners[i], corners[(i + 1) % corners.len()]);
                        boundary_normal_across_edge
                            .entry(edge)
                            .or_insert_with(Vec::new)
                            .push(*normal);
                    }
                }

                let into_part = part_pos.rotation.inverse();
                let mut boundary_faces = Vec::new();

                for ((normal, corners), boundary) in faces.iter().zip(is_boundary.iter()) {
                    if !boundary {
                        continue;
                    }

                    let edge_pseudo_normals = (0..corners.len())
                        .map(|i| {
                            let edge = edge_key(corners[i], corners[(i + 1) % corners.len()]);
                            // Halfway to the boundary face across the edge; an edge whose other
                            // face was cut away keeps this face's normal, closing the cone there.
                            let across = boundary_normal_across_edge
                                .get(&edge)
                                .and_then(|normals| {
                                    normals.iter().find(|other| **other != *normal).copied()
                                })
                                .unwrap_or(*normal);
                            into_part * (*normal + across).normalize_or(*normal)
                        })
                        .collect();

                    boundary_faces.push(CompoundFaceCone {
                        face: into_part * *normal,
                        edge_pseudo_normals,
                    });
                }

                Some(CompoundPseudoNormals { boundary_faces })
            })
            .collect();

        self.pseudo_normals = Some(pseudo_normals);
    }

    /// Returns the [`CompoundPseudoNormals`] of the part with index `i`.
    ///
    /// `None` unless this compound was given [`CompoundFlags::FIX_INTERNAL_EDGES`] and the part is
    /// polygonal -- a part with curved sides has no edge that another part could cover, so it
    /// constrains nothing.
    pub fn part_normal_constraints(&self, i: u32) -> Option<&CompoundPseudoNormals> {
        self.pseudo_normals.as_ref()?[i as usize].as_ref()
    }

    /// Creates a compound shape by decomposing a triangle mesh into convex polygons.
    ///
    /// This 2D-only method takes a [`TriMesh`] and merges adjacent triangles into larger
    /// convex polygons using the Hertel-Mehlhorn algorithm. This is useful for creating
    /// efficient collision shapes from arbitrary 2D meshes, as the resulting compound
    /// has fewer sub-shapes than using individual triangles.
    ///
    /// The algorithm works by:
    /// 1. Starting with all triangles from the input mesh
    /// 2. Merging adjacent triangles if the result is still convex
    /// 3. Creating a compound from the resulting convex polygons
    ///
    /// # Arguments
    ///
    /// * `trimesh` - A reference to the triangle mesh to decompose
    ///
    /// # Returns
    ///
    /// * `Some(Compound)` - A compound shape containing the convex polygons
    /// * `None` - If any of the created shapes has zero or near-zero area
    ///
    /// # Performance
    ///
    /// This decomposition is typically much more efficient than using the raw triangle mesh,
    /// as it reduces the number of shapes from N triangles to potentially N/2 or fewer polygons.
    ///
    /// # Example
    ///
    /// ```
    /// # #[cfg(all(feature = "dim2", feature = "f32"))] {
    /// # use parry2d::shape::{Compound, TriMesh};
    /// # use parry2d::math::Vector;
    ///
    /// // Create a simple square mesh (2 triangles)
    /// let vertices = vec![
    ///     Vector::ZERO,
    ///     Vector::new(1.0, 0.0),
    ///     Vector::new(1.0, 1.0),
    ///     Vector::new(0.0, 1.0),
    /// ];
    ///
    /// let indices = vec![
    ///     [0, 1, 2],  // First triangle
    ///     [0, 2, 3],  // Second triangle
    /// ];
    ///
    /// let trimesh = TriMesh::new(vertices, indices).unwrap();
    ///
    /// // Decompose into convex polygons
    /// if let Some(compound) = Compound::decompose_trimesh(&trimesh) {
    ///     // The two triangles should be merged into one or two convex polygons
    ///     assert!(compound.shapes().len() <= 2);
    /// }
    /// # }
    /// ```
    ///
    /// # Example: Complex Shape
    ///
    /// ```
    /// # #[cfg(all(feature = "dim2", feature = "f32"))] {
    /// # use parry2d::shape::{Compound, TriMesh};
    /// # use parry2d::math::Vector;
    ///
    /// // Create an L-shaped mesh
    /// let vertices = vec![
    ///     Vector::ZERO,
    ///     Vector::new(2.0, 0.0),
    ///     Vector::new(2.0, 1.0),
    ///     Vector::new(1.0, 1.0),
    ///     Vector::new(1.0, 2.0),
    ///     Vector::new(0.0, 2.0),
    /// ];
    ///
    /// let indices = vec![
    ///     [0, 1, 2],
    ///     [0, 2, 3],
    ///     [0, 3, 4],
    ///     [0, 4, 5],
    /// ];
    ///
    /// let trimesh = TriMesh::new(vertices, indices).unwrap();
    ///
    /// // Decompose the L-shape into convex polygons
    /// if let Some(compound) = Compound::decompose_trimesh(&trimesh) {
    ///     // The result will have fewer shapes than the original 4 triangles
    ///     assert!(compound.shapes().len() > 0);
    ///     println!("Decomposed into {} convex polygons", compound.shapes().len());
    /// }
    /// # }
    /// ```
    #[cfg(feature = "dim2")]
    pub fn decompose_trimesh(trimesh: &TriMesh) -> Option<Self> {
        let polygons = hertel_mehlhorn(trimesh.vertices(), trimesh.indices());
        let shapes: Option<Vec<_>> = polygons
            .into_iter()
            .map(|points| {
                match points.len() {
                    3 => {
                        let triangle = Triangle::new(points[0], points[1], points[2]);
                        Some(SharedShape::new(triangle))
                    }
                    _ => ConvexPolygon::from_convex_polyline(points).map(SharedShape::new),
                }
                .map(|shape| (Pose::IDENTITY, shape))
            })
            .collect();
        Some(Self::new(shapes?))
    }
}

impl Compound {
    /// Returns a slice containing all sub-shapes and their positions in this compound.
    ///
    /// Each element in the returned slice is a tuple containing:
    /// - The sub-shape's position and orientation ([`Pose`]) relative to the compound's origin
    /// - The sub-shape itself ([`SharedShape`])
    ///
    /// The order of shapes matches the order they were provided to [`Compound::new`].
    /// The index of each shape in this slice corresponds to its ID used in other operations
    /// like BVH traversal.
    ///
    /// # Returns
    ///
    /// A slice of (isometry, shape) pairs representing all sub-shapes in this compound.
    ///
    /// # Example
    ///
    /// ```
    /// # #[cfg(all(feature = "dim3", feature = "f32"))] {
    /// # use parry3d::shape::{Compound, Ball, Cuboid, SharedShape};
    /// use parry3d::math::Pose;
    /// use parry3d::math::Vector;
    ///
    /// let shapes = vec![
    ///     (Pose::translation(1.0, 0.0, 0.0), SharedShape::new(Ball::new(0.5))),
    ///     (Pose::translation(-1.0, 0.0, 0.0), SharedShape::new(Ball::new(0.3))),
    ///     (Pose::identity(), SharedShape::new(Cuboid::new(Vector::new(0.5, 0.5, 0.5)))),
    /// ];
    ///
    /// let compound = Compound::new(shapes);
    ///
    /// // Access all shapes
    /// assert_eq!(compound.shapes().len(), 3);
    ///
    /// // Inspect individual shapes
    /// for (i, (position, shape)) in compound.shapes().iter().enumerate() {
    ///     println!("Shape {} at position: {:?}", i, position.translation);
    ///
    ///     // Check if it's a ball
    ///     if let Some(ball) = shape.as_ball() {
    ///         println!("  Ball with radius: {}", ball.radius);
    ///     }
    /// }
    /// # }
    /// ```
    ///
    /// # Example: Modifying Sub-Shape Positions
    ///
    /// ```
    /// # #[cfg(all(feature = "dim3", feature = "f32"))] {
    /// # use parry3d::shape::{Compound, Ball, SharedShape};
    /// use parry3d::math::Pose;
    ///
    /// let shapes = vec![
    ///     (Pose::translation(0.0, 0.0, 0.0), SharedShape::new(Ball::new(1.0))),
    ///     (Pose::translation(2.0, 0.0, 0.0), SharedShape::new(Ball::new(1.0))),
    /// ];
    ///
    /// let compound = Compound::new(shapes);
    ///
    /// // Note: To modify positions, you need to create a new compound
    /// let mut new_shapes: Vec<_> = compound.shapes()
    ///     .iter()
    ///     .map(|(pos, shape)| {
    ///         // Shift all shapes up by 1 unit
    ///         let new_pos = pos * Pose::translation(0.0, 1.0, 0.0);
    ///         (new_pos, shape.clone())
    ///     })
    ///     .collect();
    ///
    /// let shifted_compound = Compound::new(new_shapes);
    /// assert_eq!(shifted_compound.shapes().len(), 2);
    /// # }
    /// ```
    #[inline]
    pub fn shapes(&self) -> &[(Pose, SharedShape)] {
        &self.shapes[..]
    }

    /// Returns the Axis-Aligned Bounding Box (AABB) of this compound in local space.
    ///
    /// The local AABB is the smallest axis-aligned box that contains all sub-shapes
    /// in the compound, computed in the compound's local coordinate system. This AABB
    /// is automatically computed when the compound is created and encompasses all
    /// sub-shapes at their specified positions.
    ///
    /// # Returns
    ///
    /// A reference to the [`Aabb`] representing the compound's bounding box in local space.
    ///
    /// # Use Cases
    ///
    /// - **Broad-phase collision detection**: Quick rejection tests before detailed queries
    /// - **Spatial partitioning**: Organizing compounds in larger spatial structures
    /// - **Culling**: Determining if the compound is visible or relevant to a query
    /// - **Size estimation**: Getting approximate dimensions of the compound
    ///
    /// # Example
    ///
    /// ```
    /// # #[cfg(all(feature = "dim3", feature = "f32"))] {
    /// # use parry3d::shape::{Compound, Ball, SharedShape};
    /// use parry3d::math::Pose;
    /// use parry3d::math::Vector;
    ///
    /// let shapes = vec![
    ///     (Pose::translation(-2.0, 0.0, 0.0), SharedShape::new(Ball::new(0.5))),
    ///     (Pose::translation(2.0, 0.0, 0.0), SharedShape::new(Ball::new(0.5))),
    /// ];
    ///
    /// let compound = Compound::new(shapes);
    /// let aabb = compound.local_aabb();
    ///
    /// // The AABB should contain both balls
    /// // Left ball extends from -2.5 to -1.5 on X axis
    /// // Right ball extends from 1.5 to 2.5 on X axis
    /// assert!(aabb.mins.x <= -2.5);
    /// assert!(aabb.maxs.x >= 2.5);
    ///
    /// // Check if a point is inside the AABB
    /// assert!(aabb.contains_local_point(Vector::ZERO));
    /// assert!(!aabb.contains_local_point(Vector::new(10.0, 0.0, 0.0)));
    /// # }
    /// ```
    ///
    /// # Example: Computing Compound Dimensions
    ///
    /// ```
    /// # #[cfg(all(feature = "dim3", feature = "f32"))] {
    /// # use parry3d::shape::{Compound, Cuboid, SharedShape};
    /// use parry3d::math::Pose;
    /// use parry3d::math::Vector;
    ///
    /// let shapes = vec![
    ///     (Pose::identity(), SharedShape::new(Cuboid::new(Vector::new(1.0, 1.0, 1.0)))),
    ///     (Pose::translation(3.0, 0.0, 0.0), SharedShape::new(Cuboid::new(Vector::new(0.5, 0.5, 0.5)))),
    /// ];
    ///
    /// let compound = Compound::new(shapes);
    /// let aabb = compound.local_aabb();
    ///
    /// // Calculate the total dimensions
    /// let dimensions = aabb.maxs - aabb.mins;
    /// println!("Compound dimensions: {:?}", dimensions);
    ///
    /// // The compound extends from -1.0 to 3.5 on the X axis (4.5 units total)
    /// assert!((dimensions.x - 4.5).abs() < 1e-5);
    /// # }
    /// ```
    #[inline]
    pub fn local_aabb(&self) -> &Aabb {
        &self.aabb
    }

    /// Returns the bounding sphere of this compound in local space.
    ///
    /// The bounding sphere is the smallest sphere that contains all sub-shapes in the
    /// compound. It is computed from the compound's AABB by finding the sphere that
    /// tightly encloses that box. This provides a simple, rotation-invariant bounding
    /// volume useful for certain collision detection algorithms.
    ///
    /// # Returns
    ///
    /// A [`BoundingSphere`] centered in local space that contains the entire compound.
    ///
    /// # Performance
    ///
    /// This method is very fast as it simply computes the bounding sphere from the
    /// pre-computed AABB. The bounding sphere is not cached - it's computed on each call.
    ///
    /// # Comparison with AABB
    ///
    /// - **Bounding Sphere**: Rotation-invariant, simpler intersection tests, but often looser fit
    /// - **AABB**: Tighter fit for axis-aligned objects, but must be recomputed when rotated
    ///
    /// # Example
    ///
    /// ```
    /// # #[cfg(all(feature = "dim3", feature = "f32"))] {
    /// # use parry3d::shape::{Compound, Ball, SharedShape};
    /// use parry3d::math::Pose;
    ///
    /// let shapes = vec![
    ///     (Pose::translation(-1.0, 0.0, 0.0), SharedShape::new(Ball::new(0.5))),
    ///     (Pose::translation(1.0, 0.0, 0.0), SharedShape::new(Ball::new(0.5))),
    /// ];
    ///
    /// let compound = Compound::new(shapes);
    /// let bounding_sphere = compound.local_bounding_sphere();
    ///
    /// // The bounding sphere should contain both balls
    /// println!("Center: {:?}", bounding_sphere.center());
    /// println!("Radius: {}", bounding_sphere.radius());
    ///
    /// // The center should be near the origin
    /// assert!(bounding_sphere.center().length() < 0.1);
    ///
    /// // The radius should be at least 1.5 (distance to ball edge: 1.0 + 0.5)
    /// assert!(bounding_sphere.radius() >= 1.5);
    /// # }
    /// ```
    ///
    /// # Example: Using Bounding Sphere for Quick Rejection
    ///
    /// ```
    /// # #[cfg(all(feature = "dim3", feature = "f32"))] {
    /// # use parry3d::shape::{Compound, Cuboid, SharedShape};
    /// use parry3d::math::{Pose, Vector};
    ///
    /// let shapes = vec![
    ///     (Pose::identity(), SharedShape::new(Cuboid::new(Vector::new(1.0, 1.0, 1.0)))),
    ///     (Pose::translation(2.0, 0.0, 0.0), SharedShape::new(Cuboid::new(Vector::new(0.5, 0.5, 0.5)))),
    /// ];
    ///
    /// let compound = Compound::new(shapes);
    /// let sphere = compound.local_bounding_sphere();
    ///
    /// // Quick test: is a point potentially inside the compound?
    /// let test_point = Vector::new(5.0, 5.0, 5.0);
    /// let distance_to_center = (test_point - sphere.center()).length();
    ///
    /// if distance_to_center > sphere.radius() {
    ///     println!("Vector is definitely outside the compound");
    ///     assert!(distance_to_center > sphere.radius());
    /// } else {
    ///     println!("Vector might be inside - need detailed check");
    /// }
    /// # }
    /// ```
    #[inline]
    pub fn local_bounding_sphere(&self) -> BoundingSphere {
        self.aabb.bounding_sphere()
    }

    /// Returns a slice of AABBs, one for each sub-shape in this compound.
    ///
    /// Each AABB in the returned slice corresponds to the bounding box of a sub-shape,
    /// transformed to the compound's local coordinate system. The AABBs are stored in
    /// the same order as the shapes returned by [`Compound::shapes`], so index `i` in
    /// this slice corresponds to shape `i`.
    ///
    /// These AABBs are used internally by the BVH for efficient spatial queries and
    /// collision detection. They are pre-computed during compound construction.
    ///
    /// # Returns
    ///
    /// A slice of [`Aabb`] representing the local-space bounding boxes of each sub-shape.
    ///
    /// # Use Cases
    ///
    /// - Inspecting individual sub-shape bounds without accessing the shapes themselves
    /// - Custom spatial queries or culling operations
    /// - Debugging and visualization of the compound's structure
    /// - Understanding the BVH's leaf nodes
    ///
    /// # Example
    ///
    /// ```
    /// # #[cfg(all(feature = "dim3", feature = "f32"))] {
    /// # use parry3d::shape::{Compound, Ball, Cuboid, SharedShape};
    /// use parry3d::math::Pose;
    /// use parry3d::math::Vector;
    ///
    /// let shapes = vec![
    ///     (Pose::translation(-2.0, 0.0, 0.0), SharedShape::new(Ball::new(0.5))),
    ///     (Pose::identity(), SharedShape::new(Cuboid::new(Vector::new(1.0, 1.0, 1.0)))),
    ///     (Pose::translation(3.0, 0.0, 0.0), SharedShape::new(Ball::new(0.3))),
    /// ];
    ///
    /// let compound = Compound::new(shapes);
    ///
    /// // Get AABBs for all sub-shapes
    /// let aabbs = compound.aabbs();
    /// assert_eq!(aabbs.len(), 3);
    ///
    /// // Inspect each AABB
    /// for (i, aabb) in aabbs.iter().enumerate() {
    ///     println!("Shape {} AABB:", i);
    ///     println!("  Min: {:?}", aabb.mins);
    ///     println!("  Max: {:?}", aabb.maxs);
    ///
    ///     let center = aabb.center();
    ///     let extents = aabb.half_extents();
    ///     println!("  Center: {:?}", center);
    ///     println!("  Half-extents: {:?}", extents);
    /// }
    /// # }
    /// ```
    ///
    /// # Example: Finding Sub-Shapes in a Region
    ///
    /// ```
    /// # #[cfg(all(feature = "dim3", feature = "f32"))] {
    /// # use parry3d::shape::{Compound, Ball, SharedShape};
    /// use parry3d::math::Pose;
    /// use parry3d::bounding_volume::Aabb;
    /// use parry3d::math::Vector;
    ///
    /// let shapes = vec![
    ///     (Pose::translation(-5.0, 0.0, 0.0), SharedShape::new(Ball::new(0.5))),
    ///     (Pose::translation(0.0, 0.0, 0.0), SharedShape::new(Ball::new(0.5))),
    ///     (Pose::translation(5.0, 0.0, 0.0), SharedShape::new(Ball::new(0.5))),
    /// ];
    ///
    /// let compound = Compound::new(shapes);
    ///
    /// // Define a query point
    /// let query_point = Vector::ZERO;
    ///
    /// // Find which sub-shapes might contain this point
    /// let potentially_containing: Vec<usize> = compound.aabbs()
    ///     .iter()
    ///     .enumerate()
    ///     .filter(|(_, aabb)| aabb.contains_local_point(query_point))
    ///     .map(|(i, _)| i)
    ///     .collect();
    ///
    /// // Only the middle shape (index 1) should contain the origin
    /// assert_eq!(potentially_containing.len(), 1);
    /// assert_eq!(potentially_containing[0], 1);
    /// # }
    /// ```
    #[inline]
    pub fn aabbs(&self) -> &[Aabb] {
        &self.aabbs[..]
    }

    /// Returns the Bounding Volume Hierarchy (BVH) used for efficient spatial queries.
    ///
    /// The BVH is an acceleration structure that organizes the sub-shapes hierarchically
    /// for fast collision detection and spatial queries. It enables logarithmic-time queries
    /// instead of linear searches through all sub-shapes. The BVH is automatically built
    /// when the compound is created using a binned construction strategy.
    ///
    /// # Returns
    ///
    /// A reference to the [`Bvh`] acceleration structure for this compound.
    ///
    /// # Use Cases
    ///
    /// - **Custom spatial queries**: Traverse the BVH for specialized collision detection
    /// - **Ray casting**: Efficiently find which sub-shapes intersect a ray
    /// - **AABB queries**: Find all sub-shapes intersecting a region
    /// - **Debugging**: Inspect the BVH structure and quality
    /// - **Performance analysis**: Understand query performance characteristics
    ///
    /// # Performance
    ///
    /// The BVH provides O(log n) query performance for most spatial operations, where n
    /// is the number of sub-shapes. For compounds with many shapes (100+), the BVH
    /// provides dramatic speedups compared to naive linear searches.
    ///
    /// # Example
    ///
    /// ```
    /// # #[cfg(all(feature = "dim3", feature = "f32"))] {
    /// # use parry3d::shape::{Compound, Ball, SharedShape};
    /// use parry3d::math::Pose;
    /// use parry3d::bounding_volume::Aabb;
    /// use parry3d::math::Vector;
    ///
    /// let shapes = vec![
    ///     (Pose::translation(-3.0, 0.0, 0.0), SharedShape::new(Ball::new(0.5))),
    ///     (Pose::translation(0.0, 0.0, 0.0), SharedShape::new(Ball::new(0.5))),
    ///     (Pose::translation(3.0, 0.0, 0.0), SharedShape::new(Ball::new(0.5))),
    /// ];
    ///
    /// let compound = Compound::new(shapes);
    /// let bvh = compound.bvh();
    ///
    /// // The BVH provides efficient hierarchical organization
    /// assert_eq!(bvh.leaf_count(), 3);
    /// println!("BVH root AABB: {:?}", bvh.root_aabb());
    /// # }
    /// ```
    ///
    /// # Example: Accessing BVH for Custom Queries
    ///
    /// ```
    /// # #[cfg(all(feature = "dim3", feature = "f32"))] {
    /// # use parry3d::shape::{Compound, Ball, SharedShape};
    /// use parry3d::math::Pose;
    ///
    /// let mut shapes = vec![];
    /// for i in 0..10 {
    ///     let x = i as f32 * 2.0;
    ///     shapes.push((
    ///         Pose::translation(x, 0.0, 0.0),
    ///         SharedShape::new(Ball::new(0.5))
    ///     ));
    /// }
    ///
    /// let compound = Compound::new(shapes);
    /// let bvh = compound.bvh();
    ///
    /// // The BVH organizes 10 shapes hierarchically
    /// assert_eq!(bvh.leaf_count(), 10);
    /// assert_eq!(compound.shapes().len(), 10);
    ///
    /// // Access the root AABB which bounds all shapes
    /// let root_aabb = bvh.root_aabb();
    /// println!("Root AABB spans from {:?} to {:?}", root_aabb.mins, root_aabb.maxs);
    /// # }
    /// ```
    #[inline]
    pub fn bvh(&self) -> &Bvh {
        &self.bvh
    }
}

impl CompositeShape for Compound {
    #[inline]
    fn map_part_at(
        &self,
        shape_id: u32,
        f: &mut dyn FnMut(Option<&Pose>, &dyn Shape, Option<&dyn NormalConstraints>),
    ) {
        if let Some(shape) = self.shapes.get(shape_id as usize) {
            let constraints = self
                .part_normal_constraints(shape_id)
                .map(|c| c as &dyn NormalConstraints);
            f(Some(&shape.0), &*shape.1, constraints)
        }
    }

    #[inline]
    fn bvh(&self) -> &Bvh {
        &self.bvh
    }
}

/// The identity for every value except `-0.0`, which IEEE 754 addition turns into `+0.0`. The
/// part-matching keys hash raw bit patterns, and the two zeros compare equal but hash apart.
fn without_negative_zero(coord: Real) -> Real {
    coord + 0.0
}

/// The edges of a closed outline, as consecutive corner pairs wrapping back to the start.
#[cfg(feature = "dim2")]
fn outline_edges(outline: &[Vector]) -> impl Iterator<Item = (Vector, Vector)> + '_ {
    outline
        .iter()
        .enumerate()
        .map(move |(i, corner)| (*corner, outline[(i + 1) % outline.len()]))
}

/// Twice the area enclosed by a closed outline, negative when it is wound clockwise.
#[cfg(feature = "dim2")]
fn signed_area(outline: &[Vector]) -> Real {
    outline_edges(outline)
        .map(|(a, b)| a.x * b.y - a.y * b.x)
        .sum()
}

impl TypedCompositeShape for Compound {
    type PartShape = dyn Shape;
    type PartNormalConstraints = CompoundPseudoNormals;

    #[inline(always)]
    fn map_typed_part_at<T>(
        &self,
        i: u32,
        mut f: impl FnMut(Option<&Pose>, &Self::PartShape, Option<&Self::PartNormalConstraints>) -> T,
    ) -> Option<T> {
        let (part_pos, part) = &self.shapes[i as usize];
        Some(f(Some(part_pos), &**part, self.part_normal_constraints(i)))
    }

    #[inline(always)]
    fn map_untyped_part_at<T>(
        &self,
        i: u32,
        mut f: impl FnMut(Option<&Pose>, &Self::PartShape, Option<&dyn NormalConstraints>) -> T,
    ) -> Option<T> {
        let (part_pos, part) = &self.shapes[i as usize];
        let constraints = self
            .part_normal_constraints(i)
            .map(|c| c as &dyn NormalConstraints);
        Some(f(Some(part_pos), &**part, constraints))
    }
}

#[cfg(test)]
#[cfg(all(feature = "dim2", feature = "alloc", feature = "std"))]
mod test {
    use super::{Compound, CompoundFlags};
    use crate::math::{Pose, Real, Vector};
    use crate::query::{ContactManifold, DefaultQueryDispatcher, PersistentQueryDispatcher};
    use crate::shape::{ConvexPolygon, Cuboid, SharedShape, Triangle};
    use alloc::{vec, vec::Vec};

    /// A flat run that turns into a ramp, decomposed the way a convex decomposition would leave it:
    /// a big piece below and a wedge above, meeting along a cut that runs right from the ramp foot
    /// at the origin. The ramp face and the cut leave the wedge with a tip pointing back along the
    /// flat run, which is what a body slides into.
    fn ramp_foot() -> Compound {
        let part = |points: &[[Real; 2]]| {
            let points = points.iter().map(|p| Vector::new(p[0], p[1])).collect();
            (
                Pose::IDENTITY,
                SharedShape::new(
                    ConvexPolygon::from_convex_polyline_unmodified(points).expect("convex part"),
                ),
            )
        };

        Compound::new(vec![
            part(&[
                [-3.825, -0.015],
                [-11.913, -6.486],
                [5.166, -6.486],
                [5.165, -0.015],
                [0.0, 0.0],
            ]),
            part(&[[0.0, 0.0], [5.165, -0.015], [3.432, 1.912], [2.135, 1.891]]),
        ])
    }

    /// The normals the wedge reports against a box whose right edge sits `overshoot` past the tip.
    fn normals_against_the_wedge(compound: &Compound, overshoot: Real) -> Vec<Vector> {
        let slider = Cuboid::new(Vector::new(0.42, 0.20));
        let pos12 = Pose::from_parts(Vector::new(overshoot - 0.42, 0.198), Default::default());

        let mut manifolds: Vec<ContactManifold<(), ()>> = Vec::new();
        let mut workspace = None;
        DefaultQueryDispatcher
            .contact_manifolds(
                &pos12,
                compound,
                &slider,
                0.02,
                &mut manifolds,
                &mut workspace,
            )
            .expect("dispatch");

        manifolds
            .iter()
            .filter(|manifold| manifold.subshape1 == 1 && !manifold.points.is_empty())
            .map(|manifold| manifold.local_n1)
            .collect()
    }

    #[test]
    fn cut_edges_are_dropped_and_close_the_corner_they_end_at() {
        let mut compound = ramp_foot();
        compound.set_flags(CompoundFlags::FIX_INTERNAL_EDGES);

        // Each part loses exactly the cut it shares with the other.
        assert_eq!(
            compound
                .part_normal_constraints(0)
                .unwrap()
                .boundary_edges
                .len(),
            4
        );
        let wedge = compound.part_normal_constraints(1).unwrap();
        assert_eq!(wedge.boundary_edges.len(), 3);

        // The ramp face runs into the cut at the tip, so its cone closes onto itself there.
        let ramp = wedge
            .boundary_edges
            .iter()
            .find(|cone| cone.face.x < -0.5)
            .expect("ramp face");
        assert_eq!(ramp.counter_clockwise_limit, ramp.face);
    }

    /// Nothing guarantees which way a part was wound, and a clockwise one would have every normal
    /// pointing into itself.
    #[test]
    fn a_clockwise_part_still_faces_outwards() {
        let clockwise = Triangle::new(
            Vector::new(0.0, 0.0),
            Vector::new(0.0, 1.0),
            Vector::new(1.0, 0.0),
        );
        let mut compound = Compound::new(vec![(Pose::IDENTITY, SharedShape::new(clockwise))]);
        compound.set_flags(CompoundFlags::FIX_INTERNAL_EDGES);

        let cones = compound.part_normal_constraints(0).expect("polygonal part");
        assert_eq!(cones.boundary_edges.len(), 3);

        // Every face points away from the triangle's interior, which lies towards the origin.
        let interior = Vector::new(0.25, 0.25);
        for cone in &cones.boundary_edges {
            assert!(
                cone.face.dot(interior) < 0.0 || cone.face.dot(Vector::new(1.0, 1.0)) > 0.0,
                "inward face {:?}",
                cone.face
            );
        }

        // The hypotenuse is the only face pointing up and to the right.
        let hypotenuse = cones
            .boundary_edges
            .iter()
            .find(|cone| cone.face.x > 0.5 && cone.face.y > 0.5)
            .expect("outward hypotenuse");
        assert!(hypotenuse.face.dot(Vector::new(1.0, 1.0).normalize()) > 0.99);
    }

    #[test]
    fn a_body_running_into_the_tip_reports_no_contact_with_the_cut() {
        let plain = ramp_foot();
        let mut fixed = ramp_foot();
        fixed.set_flags(CompoundFlags::FIX_INTERNAL_EDGES);

        // Against the bare compound the wedge answers with its tip axis: a wall across the flat run.
        let ghost = normals_against_the_wedge(&plain, 0.0);
        assert_eq!(ghost.len(), 1);
        assert!(
            ghost[0].x < -0.99,
            "expected a wall normal, got {:?}",
            ghost[0]
        );

        // The cut cannot be touched, so the wedge reports nothing until the ramp itself is reached.
        assert!(normals_against_the_wedge(&fixed, 0.0).is_empty());
    }

    #[test]
    fn queries_identify_the_part_they_hit() {
        use crate::query::{PointQuery, Ray, RayCast};
        use crate::shape::FeatureId;

        let compound = ramp_foot();

        // Straight down onto the flat run: the big bottom piece, part 0.
        let onto_the_flat = Ray::new(Vector::new(-2.0, 5.0), Vector::new(0.0, -1.0));
        let hit = compound
            .cast_local_ray_and_get_normal(&onto_the_flat, Real::MAX, true)
            .expect("hits the flat");
        assert_eq!(hit.feature, FeatureId::Face(0));

        // Straight down onto the ramp: the wedge, part 1.
        let onto_the_ramp = Ray::new(Vector::new(3.0, 5.0), Vector::new(0.0, -1.0));
        let hit = compound
            .cast_local_ray_and_get_normal(&onto_the_ramp, Real::MAX, true)
            .expect("hits the ramp");
        assert_eq!(hit.feature, FeatureId::Face(1));

        let (_, feature) = compound.project_local_point_and_get_feature(Vector::new(3.0, 3.0));
        assert_eq!(feature, FeatureId::Face(1));
    }

    #[test]
    fn contacts_with_the_ramp_itself_are_left_alone() {
        let plain = ramp_foot();
        let mut fixed = ramp_foot();
        fixed.set_flags(CompoundFlags::FIX_INTERNAL_EDGES);

        let expected = normals_against_the_wedge(&plain, 0.01);
        assert_eq!(expected.len(), 1);
        assert!(
            expected[0].y > 0.7,
            "expected the ramp face, got {:?}",
            expected[0]
        );
        assert_eq!(normals_against_the_wedge(&fixed, 0.01), expected);
    }
}

#[cfg(test)]
#[cfg(all(feature = "dim3", feature = "alloc", feature = "std"))]
mod dim3_test {
    use super::{Compound, CompoundFlags};
    use crate::math::{Pose, Real, Vector};
    use crate::shape::SharedShape;
    use alloc::{vec, vec::Vec};

    /// A flat run turning into a ramp, split along the diagonal a convex decomposition would cut:
    /// a slab and a wedge, extruded, sharing one quad exactly. Built through `convex_hull`, the
    /// way the Godot plugin builds every `ConvexPolygonShape3D`.
    fn ramp_foot() -> Compound {
        let prism = |outline: &[[Real; 2]]| {
            let corners: Vec<Vector> = outline
                .iter()
                .flat_map(|p| [Vector::new(p[0], p[1], -3.0), Vector::new(p[0], p[1], 3.0)])
                .collect();
            (
                Pose::IDENTITY,
                SharedShape::convex_hull(&corners).expect("convex prism"),
            )
        };

        Compound::new(vec![
            prism(&[[-6.0, 0.0], [5.0, 0.0], [11.0, -6.0], [-6.0, -6.0]]),
            prism(&[[5.0, 0.0], [9.0, 3.0], [11.0, 3.0], [11.0, -6.0]]),
        ])
    }

    #[test]
    fn cut_faces_are_dropped_from_both_parts() {
        let mut compound = ramp_foot();
        compound.set_flags(CompoundFlags::FIX_INTERNAL_EDGES);

        // Each prism has 6 faces and loses the shared quad -- which only holds if `convex_hull`
        // produced the cut with the same corners on both sides, the premise everything rests on.
        for part in 0..2 {
            let cones = compound.part_normal_constraints(part).expect("polyhedron");
            assert_eq!(cones.boundary_faces.len(), 5, "part {part}");
        }
    }

    /// A body pressed into the wedge tip reports a normal along its own face, not a surface of
    /// the union; the cones pull it onto the one face that is.
    #[test]
    fn a_tip_normal_is_clamped_onto_the_ramp_face() {
        use crate::query::details::NormalConstraints;

        let mut compound = ramp_foot();
        compound.set_flags(CompoundFlags::FIX_INTERNAL_EDGES);
        let wedge = compound.part_normal_constraints(1).expect("polyhedron");

        let ramp = Vector::new(-0.6, 0.8, 0.0);
        let clamped = wedge
            .project_local_normal(Vector::new(-1.0, 0.0, 0.0))
            .expect("tip direction faces the outline");
        assert!((clamped - ramp).length() < 1.0e-4, "clamped to {clamped:?}");

        // A direction the outline genuinely faces is left alone.
        let kept = wedge.project_local_normal(ramp).expect("faces the outline");
        assert!((kept - ramp).length() < 1.0e-4, "moved to {kept:?}");
    }
}
