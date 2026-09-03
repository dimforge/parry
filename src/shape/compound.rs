//!
//! Shape composed from the union of primitives.
//!

use crate::bounding_volume::{Aabb, BoundingSphere, BoundingVolume};
use crate::math::Pose;
use crate::math::{Real, Vector};
use crate::partitioning::{Bvh, BvhBuildStrategy};
use crate::query::details::NormalConstraints;
#[cfg(feature = "dim2")]
use crate::shape::compound_pseudo_normals::turns_clockwise;
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
    /// The tolerance [`Compound::set_flags`] welds part corners with, in ULPs.
    ///
    /// Two spellings of the same corner drift by about one ULP per rounding, and a corner reaches
    /// the surface through a rotation and a translation, so a handful of ULPs covers the seams a
    /// shared grid produces while staying far below any real feature.
    pub const DEFAULT_WELD_TOLERANCE: Real = 4.0;

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

    /// Builds a new compound shape from a collection of sub-shapes, with the given
    /// [`CompoundFlags`] already applied.
    ///
    /// `weld_tolerance` is [`Compound::set_flags`]'s, and `None` selects
    /// [`Compound::DEFAULT_WELD_TOLERANCE`].
    ///
    /// # Example
    ///
    /// ```rust
    /// # #[cfg(all(feature = "dim2", feature = "f32"))] {
    /// use parry2d::math::{Pose, Vector};
    /// use parry2d::shape::{Compound, CompoundFlags, Cuboid, SharedShape};
    ///
    /// // Two boxes meeting along one edge, as a decomposition or a voxel grid leaves them.
    /// let box_at = |x: f32| {
    ///     (
    ///         Pose::from_parts(Vector::new(x, 0.0), Default::default()),
    ///         SharedShape::new(Cuboid::new(Vector::new(0.5, 0.5))),
    ///     )
    /// };
    /// let compound = Compound::with_flags(
    ///     vec![box_at(0.0), box_at(1.0)],
    ///     CompoundFlags::FIX_INTERNAL_EDGES,
    ///     None,
    /// );
    ///
    /// // The edge they share is interior to the union, so neither part reports it.
    /// assert_eq!(compound.part_normal_constraints(0).unwrap().boundary_edges.len(), 3);
    /// # }
    /// ```
    pub fn with_flags(
        shapes: Vec<(Pose, SharedShape)>,
        flags: CompoundFlags,
        weld_tolerance: Option<Real>,
    ) -> Compound {
        let mut compound = Self::new(shapes);
        compound.set_flags(flags, weld_tolerance);
        compound
    }

    /// Sets the [`CompoundFlags`], computing or discarding the compound's optional associated data.
    ///
    /// `weld_tolerance` says how far apart two part corners may be and still be treated as the same
    /// point, in ULPs; `None` selects [`Compound::DEFAULT_WELD_TOLERANCE`].
    ///
    /// Parts of a decomposition meet on corners each computes for itself, and two spellings of the
    /// same corner need not round to the same bits: a row of boxes on a shared grid, or any part
    /// placed through a rotation, routinely leaves its seams a few ULPs apart. Corners that far
    /// apart are welded so the edge or face between them is still recognized as shared.
    ///
    /// The tolerance is relative to each corner's own distance from the origin, not the whole
    /// compound's, so detail near the origin keeps a tolerance of its own magnitude however far the
    /// rest of the compound reaches. Zero leaves only corners at the very same coordinates welded.
    ///
    /// Only [`CompoundFlags::FIX_INTERNAL_EDGES`] reads the tolerance; without it nothing is
    /// computed and the value is ignored.
    pub fn set_flags(&mut self, flags: CompoundFlags, weld_tolerance: Option<Real>) {
        self.flags = flags;

        if flags.contains(CompoundFlags::FIX_INTERNAL_EDGES) {
            self.compute_pseudo_normals(weld_tolerance.unwrap_or(Self::DEFAULT_WELD_TOLERANCE));
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
    /// An edge two parts share is a cut, so it is dropped, and each surviving edge opens towards
    /// the edge next to it along the union's outline, which a cut may have left in another part.
    /// That closes the corner where a surface runs into a cut -- the ledge a body would otherwise
    /// catch on -- without closing one the cut merely passes through.
    #[cfg(feature = "dim2")]
    fn compute_pseudo_normals(&mut self, weld_tolerance: Real) {
        let outlines = self.part_outlines();

        // Weld every part's corners into one vertex numbering, so two parts meeting along an edge
        // key it on the same pair of vertices even though each rounded its own corners.
        let corners: Vec<Vector> = outlines.iter().flatten().flatten().copied().collect();
        let welded = weld_corners(&corners, weld_tolerance);

        let mut vertex_ids: Vec<Option<Vec<u32>>> = Vec::with_capacity(outlines.len());
        let mut next = 0;
        for outline in &outlines {
            match outline {
                Some(outline) => {
                    vertex_ids.push(Some(welded[next..next + outline.len()].to_vec()));
                    next += outline.len();
                }
                None => vertex_ids.push(None),
            }
        }

        let edge_key = |a: u32, b: u32| if a <= b { (a, b) } else { (b, a) };

        let mut parts_sharing_edge = HashMap::default();
        for ids in vertex_ids.iter().flatten() {
            for (a, b) in outline_edge_ids(ids) {
                *parts_sharing_edge.entry(edge_key(a, b)).or_insert(0u32) += 1;
            }
        }

        // Per part, the outward normal of every outline edge, or `None` where a sibling part covers
        // it: that edge is a cut, interior to the union.
        let boundary_normals: Vec<Option<Vec<Option<Vector>>>> = outlines
            .iter()
            .zip(vertex_ids.iter())
            .map(|(outline, ids)| {
                let (outline, ids) = (outline.as_ref()?, ids.as_ref()?);

                Some(
                    outline_edges(outline)
                        .zip(outline_edge_ids(ids))
                        .map(|((a, b), (id_a, id_b))| {
                            let shared = parts_sharing_edge
                                .get(&edge_key(id_a, id_b))
                                .is_some_and(|parts| *parts > 1);

                            (!shared)
                                .then(|| crate::utils::ccw_face_normal([a, b]))
                                .flatten()
                        })
                        .collect(),
                )
            })
            .collect();

        // The union's boundary edge arriving at, and the one leaving, each corner. Gathered across
        // every part, so a corner a cut splits still finds the neighbour that continues the
        // outline into the sibling part. `None` marks a corner more than one edge claims, where no
        // single neighbour can be picked.
        let mut arriving = HashMap::default();
        let mut leaving = HashMap::default();
        for (ids, normals) in vertex_ids.iter().zip(boundary_normals.iter()) {
            let (Some(ids), Some(normals)) = (ids.as_ref(), normals.as_ref()) else {
                continue;
            };

            for ((a, b), normal) in outline_edge_ids(ids).zip(normals.iter()) {
                let Some(normal) = *normal else {
                    continue;
                };

                let _ = arriving
                    .entry(b)
                    .and_modify(|slot| *slot = None)
                    .or_insert(Some(normal));
                let _ = leaving
                    .entry(a)
                    .and_modify(|slot| *slot = None)
                    .or_insert(Some(normal));
            }
        }

        let pseudo_normals = self
            .shapes
            .iter()
            .zip(vertex_ids.iter())
            .zip(boundary_normals.iter())
            .map(|(((part_pos, _), ids), normals)| {
                let (ids, normals) = (ids.as_ref()?, normals.as_ref()?);

                let into_part = part_pos.rotation.inverse();
                let mut boundary_edges = Vec::new();

                for ((a, b), face) in outline_edge_ids(ids).zip(normals.iter()) {
                    let Some(face) = *face else {
                        continue;
                    };

                    // Halfway, so the two edges at a corner split its range and neither cone grows
                    // past a half turn, which the containment test could not represent.
                    let halfway_to = |neighbour: Option<Vector>| match neighbour {
                        Some(neighbour) => (face + neighbour).normalize_or(face),
                        None => face,
                    };

                    // A concave corner has no outward range for its two edges to share, so the cone
                    // stops on this edge's own normal, as it does where the surface runs into a cut.
                    // Walking counter-clockwise turns the normals the same way, so the edge arriving
                    // at a corner lies clockwise of the one leaving it.
                    let clockwise = arriving
                        .get(&a)
                        .copied()
                        .flatten()
                        .filter(|previous| !turns_clockwise(*previous, face));
                    let counter_clockwise = leaving
                        .get(&b)
                        .copied()
                        .flatten()
                        .filter(|next| !turns_clockwise(face, *next));

                    // Cones are stated in the part's frame, which is where the narrow phase applies
                    // them.
                    boundary_edges.push(CompoundEdgeCone {
                        face: into_part * face,
                        clockwise_limit: into_part * halfway_to(clockwise),
                        counter_clockwise_limit: into_part * halfway_to(counter_clockwise),
                    });
                }

                Some(CompoundPseudoNormals { boundary_edges })
            })
            .collect();

        self.pseudo_normals = Some(pseudo_normals);
    }

    /// The parts' faces in the compound's frame, each an outward normal with its vertex loop wound
    /// counter-clockwise around it. `None` for a part with no flat faces, which therefore shares
    /// none.
    #[cfg(feature = "dim3")]
    #[allow(clippy::type_complexity)]
    fn part_faces(&self) -> Vec<Option<Vec<(Vector, Vec<Vector>)>>> {
        self.shapes
            .iter()
            .map(|(part_pos, part)| {
                let faces = if let Some(polyhedron) = part.as_convex_polyhedron() {
                    let points = polyhedron.points();
                    let corner_indices = polyhedron.vertices_adj_to_face();

                    polyhedron
                        .faces()
                        .iter()
                        .map(|face| {
                            let first = face.first_vertex_or_edge as usize;
                            let corners = corner_indices
                                [first..first + face.num_vertices_or_edges as usize]
                                .iter()
                                .map(|index| points[*index as usize])
                                .collect();
                            (face.normal, corners)
                        })
                        .collect()
                } else {
                    // A cuboid is a polyhedron the type system does not call one, and a row of them
                    // is how a voxel grid reaches a compound.
                    cuboid_faces(part.as_cuboid()?.half_extents)
                };

                Some(
                    faces
                        .into_iter()
                        .map(|(normal, corners): (Vector, Vec<Vector>)| {
                            (
                                part_pos.rotation * normal,
                                corners
                                    .into_iter()
                                    .map(|corner| part_pos.transform_point(corner))
                                    .collect(),
                            )
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
    /// towards the boundary face across it, which a cut may have left in another part. That closes
    /// the corner where a surface runs into a cut -- the ledge a body would otherwise catch on --
    /// without closing one the cut merely passes through.
    #[cfg(feature = "dim3")]
    fn compute_pseudo_normals(&mut self, weld_tolerance: Real) {
        let part_faces = self.part_faces();

        // Weld every part's corners into one vertex numbering, so two parts meeting along a face
        // key it on the same set of vertices even though each rounded its own corners.
        let corners: Vec<Vector> = part_faces
            .iter()
            .flatten()
            .flatten()
            .flat_map(|(_, corners)| corners.iter().copied())
            .collect();
        let welded = weld_corners(&corners, weld_tolerance);

        let mut vertex_ids: Vec<Option<Vec<Vec<u32>>>> = Vec::with_capacity(part_faces.len());
        let mut next = 0;
        for faces in &part_faces {
            match faces {
                Some(faces) => {
                    let mut per_face = Vec::with_capacity(faces.len());
                    for (_, corners) in faces {
                        per_face.push(welded[next..next + corners.len()].to_vec());
                        next += corners.len();
                    }
                    vertex_ids.push(Some(per_face));
                }
                None => vertex_ids.push(None),
            }
        }

        // A shared face has the same vertices on both sides, in some rotation and possibly reversed,
        // so the key ignores their order.
        let face_key = |ids: &[u32]| {
            let mut key = ids.to_vec();
            key.sort_unstable();
            key
        };
        let edge_key = |a: u32, b: u32| if a <= b { (a, b) } else { (b, a) };

        let mut parts_sharing_face = HashMap::default();
        for ids in vertex_ids.iter().flatten() {
            for face in ids {
                *parts_sharing_face.entry(face_key(face)).or_insert(0u32) += 1;
            }
        }

        // Per part, whether each face lies on the union's outline rather than being a cut two
        // parts share.
        let is_boundary: Vec<Option<Vec<bool>>> = vertex_ids
            .iter()
            .map(|ids| {
                Some(
                    ids.as_ref()?
                        .iter()
                        .map(|face| {
                            parts_sharing_face
                                .get(&face_key(face))
                                .is_none_or(|parts| *parts <= 1)
                        })
                        .collect(),
                )
            })
            .collect();

        // The boundary faces meeting along each edge. Gathered across every part, so an edge a cut
        // splits still finds the face continuing the surface into the sibling part.
        let mut faces_along_edge = HashMap::default();
        for ((faces, ids), boundary) in part_faces
            .iter()
            .zip(vertex_ids.iter())
            .zip(is_boundary.iter())
        {
            let (Some(faces), Some(ids), Some(boundary)) =
                (faces.as_ref(), ids.as_ref(), boundary.as_ref())
            else {
                continue;
            };

            for (((normal, _), face), boundary) in faces.iter().zip(ids.iter()).zip(boundary.iter())
            {
                if !boundary {
                    continue;
                }

                for i in 0..face.len() {
                    faces_along_edge
                        .entry(edge_key(face[i], face[(i + 1) % face.len()]))
                        .or_insert_with(Vec::new)
                        .push(*normal);
                }
            }
        }

        let pseudo_normals = self
            .shapes
            .iter()
            .zip(part_faces.iter())
            .zip(vertex_ids.iter())
            .zip(is_boundary.iter())
            .map(|((((part_pos, _), faces), ids), boundary)| {
                let (faces, ids, boundary) = (faces.as_ref()?, ids.as_ref()?, boundary.as_ref()?);

                let into_part = part_pos.rotation.inverse();
                let mut boundary_faces = Vec::new();

                for (((normal, corners), face), on_boundary) in
                    faces.iter().zip(ids.iter()).zip(boundary.iter())
                {
                    if !on_boundary {
                        continue;
                    }

                    let edge_pseudo_normals = (0..corners.len())
                        .map(|i| {
                            let next = (i + 1) % corners.len();
                            let (from, to) = (corners[i], corners[next]);
                            // Halfway to the boundary face across the edge, which may belong to
                            // another part where a cut splits the surface along it. A concave edge
                            // has no outward range for its two faces to share, so the cone stops on
                            // this face's own normal, as it does at an edge whose other face was
                            // cut away. Faces are wound counter-clockwise around their normal, so
                            // the edge is convex exactly when the two normals turn that way too.
                            let across = faces_along_edge
                                .get(&edge_key(face[i], face[next]))
                                .and_then(|normals| {
                                    normals.iter().find(|other| **other != *normal).copied()
                                })
                                .filter(|across| normal.cross(*across).dot(to - from) >= 0.0)
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

/// The six faces of a cuboid in its own frame, each an outward normal with its corner loop wound
/// counter-clockwise around it.
#[cfg(feature = "dim3")]
#[allow(clippy::type_complexity)]
fn cuboid_faces(half_extents: Vector) -> Vec<(Vector, Vec<Vector>)> {
    let (x, y, z) = (half_extents.x, half_extents.y, half_extents.z);
    let corner = |sx: Real, sy: Real, sz: Real| Vector::new(sx * x, sy * y, sz * z);

    alloc::vec![
        (
            Vector::X,
            alloc::vec![
                corner(1.0, -1.0, -1.0),
                corner(1.0, 1.0, -1.0),
                corner(1.0, 1.0, 1.0),
                corner(1.0, -1.0, 1.0),
            ]
        ),
        (
            -Vector::X,
            alloc::vec![
                corner(-1.0, -1.0, -1.0),
                corner(-1.0, -1.0, 1.0),
                corner(-1.0, 1.0, 1.0),
                corner(-1.0, 1.0, -1.0),
            ]
        ),
        (
            Vector::Y,
            alloc::vec![
                corner(-1.0, 1.0, -1.0),
                corner(-1.0, 1.0, 1.0),
                corner(1.0, 1.0, 1.0),
                corner(1.0, 1.0, -1.0),
            ]
        ),
        (
            -Vector::Y,
            alloc::vec![
                corner(-1.0, -1.0, -1.0),
                corner(1.0, -1.0, -1.0),
                corner(1.0, -1.0, 1.0),
                corner(-1.0, -1.0, 1.0),
            ]
        ),
        (
            Vector::Z,
            alloc::vec![
                corner(-1.0, -1.0, 1.0),
                corner(1.0, -1.0, 1.0),
                corner(1.0, 1.0, 1.0),
                corner(-1.0, 1.0, 1.0),
            ]
        ),
        (
            -Vector::Z,
            alloc::vec![
                corner(-1.0, -1.0, -1.0),
                corner(-1.0, 1.0, -1.0),
                corner(1.0, 1.0, -1.0),
                corner(1.0, -1.0, -1.0),
            ]
        ),
    ]
}

/// Assigns every corner the index of the vertex it welds onto, so parts that meet are keyed on
/// the same point even when each computed it for itself.
///
/// Two corners weld when no coordinate differs by more than `tolerance` ULPs, measured at the
/// scale of whichever of the two lies further from the origin. That scale has to come from the
/// point rather than the coordinate: a corner reaching the outline through a rotation can land
/// near an axis, leaving a coordinate whose own magnitude is pure cancellation and says nothing
/// about the error in it, while its siblings still carry the magnitude the rounding happened at.
///
/// A tolerance of zero leaves only corners at the very same coordinates welded, matching what a
/// shared vertex buffer produces and nothing else. The comparison stays numeric there rather than
/// dropping to raw bits, so the two signed zeros still meet.
fn weld_corners(corners: &[Vector], tolerance: Real) -> Vec<u32> {
    let radius = |corner: &Vector| tolerance * Real::EPSILON * corner.abs().max_element();

    // Union-find over the corners, so a chain of near-coincident corners collapses to one vertex
    // whichever order the pairs turn up in.
    let mut parent: Vec<u32> = (0..corners.len() as u32).collect();
    fn find(parent: &mut [u32], mut i: u32) -> u32 {
        while parent[i as usize] != i {
            parent[i as usize] = parent[parent[i as usize] as usize];
            i = parent[i as usize];
        }
        i
    }

    let leaves: Vec<Aabb> = corners
        .iter()
        .map(|corner| Aabb::new(*corner, *corner))
        .collect();
    let bvh = Bvh::from_leaves(BvhBuildStrategy::Binned, &leaves);

    for (i, corner) in corners.iter().enumerate() {
        // Querying every corner with its own radius finds every pair that should weld: the test is
        // symmetric and scaled by the larger of the two, so the larger corner's own query reaches
        // the smaller one even when the smaller one's radius would fall short.
        let reach = radius(corner);
        let query = Aabb::new(
            *corner - Vector::splat(reach),
            *corner + Vector::splat(reach),
        );

        for j in bvh.intersect_aabb(&query) {
            let other = corners[j as usize];
            let scale = corner.abs().max_element().max(other.abs().max_element());

            if (*corner - other).abs().max_element() <= tolerance * Real::EPSILON * scale {
                let (a, b) = (find(&mut parent, i as u32), find(&mut parent, j));
                parent[a as usize] = b;
            }
        }
    }

    // Number the roots so the ids are dense.
    let mut ids = HashMap::default();
    (0..corners.len() as u32)
        .map(|i| {
            let root = find(&mut parent, i);
            let next = ids.len() as u32;
            *ids.entry(root).or_insert(next)
        })
        .collect()
}

/// The edges of a closed outline, as consecutive vertex-index pairs wrapping back to the start.
#[cfg(feature = "dim2")]
fn outline_edge_ids(ids: &[u32]) -> impl Iterator<Item = (u32, u32)> + '_ {
    ids.iter()
        .enumerate()
        .map(move |(i, id)| (*id, ids[(i + 1) % ids.len()]))
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
    use crate::math::{Pose, Real, Rotation, Vector};
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
        compound.set_flags(CompoundFlags::FIX_INTERNAL_EDGES, None);

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
        compound.set_flags(CompoundFlags::FIX_INTERNAL_EDGES, None);

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
        fixed.set_flags(CompoundFlags::FIX_INTERNAL_EDGES, None);

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

    /// An L wound counter-clockwise, split along the diagonal a convex decomposition would cut
    /// from its reflex corner. The cut lands on two corners of the union: the convex one at the
    /// origin, and the reflex one at `(1, 1)`.
    fn cut_ell() -> Compound {
        let part = |points: &[[Real; 2]]| {
            let points = points.iter().map(|p| Vector::new(p[0], p[1])).collect();
            (
                Pose::IDENTITY,
                SharedShape::new(
                    ConvexPolygon::from_convex_polyline_unmodified(points).expect("convex part"),
                ),
            )
        };

        let mut compound = Compound::new(vec![
            part(&[[0.0, 0.0], [2.0, 0.0], [2.0, 1.0], [1.0, 1.0]]),
            part(&[[0.0, 0.0], [1.0, 1.0], [1.0, 2.0], [0.0, 2.0]]),
        ]);
        compound.set_flags(CompoundFlags::FIX_INTERNAL_EDGES, None);
        compound
    }

    /// The cone of the boundary edge of `part` whose outward normal is `face`.
    fn cone_facing(compound: &Compound, part: u32, face: Vector) -> super::CompoundEdgeCone {
        *compound
            .part_normal_constraints(part)
            .expect("polygonal part")
            .boundary_edges
            .iter()
            .find(|cone| (cone.face - face).length() < 1.0e-5)
            .unwrap_or_else(|| panic!("no edge facing {face:?}"))
    }

    /// A corner a cut splits is still a corner of the union, so the two edges meeting there open
    /// towards each other even though they belong to different parts.
    #[test]
    fn a_convex_corner_split_by_a_cut_opens_across_the_parts() {
        let compound = cut_ell();

        // Each part keeps its three outline edges and loses the diagonal it shares.
        for part in 0..2 {
            assert_eq!(
                compound
                    .part_normal_constraints(part)
                    .unwrap()
                    .boundary_edges
                    .len(),
                3,
                "part {part}"
            );
        }

        // The union turns 90 degrees at the origin, between the bottom edge of one part and the
        // left edge of the other. Both cones reach the bisector, splitting the corner between them.
        let bisector = Vector::new(-1.0, -1.0).normalize();
        let bottom = cone_facing(&compound, 0, -Vector::Y);
        let left = cone_facing(&compound, 1, -Vector::X);

        assert!(
            (bottom.clockwise_limit - bisector).length() < 1.0e-5,
            "bottom stopped at {:?}",
            bottom.clockwise_limit
        );
        assert!(
            (left.counter_clockwise_limit - bisector).length() < 1.0e-5,
            "left stopped at {:?}",
            left.counter_clockwise_limit
        );
    }

    /// A concave corner has no outward range for its edges to share, so each keeps its own normal
    /// -- the same thing that happens where a surface runs into a cut.
    #[test]
    fn a_concave_corner_split_by_a_cut_stays_closed() {
        let compound = cut_ell();

        // The reflex corner of the L, at `(1, 1)`: the upwards edge of one part meets the
        // rightwards edge of the other.
        let upwards = cone_facing(&compound, 0, Vector::Y);
        let rightwards = cone_facing(&compound, 1, Vector::X);

        assert_eq!(upwards.counter_clockwise_limit, upwards.face);
        assert_eq!(rightwards.clockwise_limit, rightwards.face);
    }

    /// A row of boxes on a shared grid, the way a voxel body reaches a compound, optionally
    /// rotated as a whole. Every part computes its own corners, so the seams are only found if
    /// corners that round differently still weld.
    fn voxel_row(n: u32, size: Real, origin: Real, angle: Real) -> Compound {
        let half = size / 2.0;
        let rotation = Rotation::from_angle(angle);
        let parts = (0..n)
            .map(|i| {
                let center = Vector::new(origin + i as Real * size, 0.0);
                (
                    Pose::from_parts(rotation * center, rotation),
                    SharedShape::new(Cuboid::new(Vector::new(half, half))),
                )
            })
            .collect();

        let mut compound = Compound::new(parts);
        compound.set_flags(CompoundFlags::FIX_INTERNAL_EDGES, None);
        compound
    }

    /// How many boundary edges each part of a row kept. A welded row leaves the two ends with
    /// three and every box between them with two.
    fn boundary_edge_counts(compound: &Compound, n: u32) -> Vec<usize> {
        (0..n)
            .map(|i| {
                compound
                    .part_normal_constraints(i)
                    .map_or(99, |cones| cones.boundary_edges.len())
            })
            .collect()
    }

    /// Parts of a grid meet on corners each computes for itself, and the two spellings need not
    /// round alike. Every row here leaves at least one seam a few ULPs wide.
    #[test]
    fn a_voxel_row_finds_seams_its_parts_rounded_differently() {
        for &(size, origin, angle) in &[
            (1.0 as Real, 0.0 as Real, 0.1 as Real),
            (1.0, 0.0, 0.5236),
            (1.0, 0.0, 0.7854),
            (0.1, 0.0, 0.0),
            (0.7, 0.0, 0.0),
            (1.0, 0.05, 0.0),
            (0.1, 1000.0, 0.0),
            (1.0, 1000.0, 0.7854),
        ] {
            assert_eq!(
                boundary_edge_counts(&voxel_row(5, size, origin, angle), 5),
                alloc::vec![3, 2, 2, 2, 3],
                "size {size} origin {origin} angle {angle}"
            );
        }
    }

    /// A zero tolerance is the exact comparison, which only recognizes a seam whose two spellings
    /// round to the very same bits.
    #[test]
    fn a_zero_tolerance_welds_only_bit_identical_corners() {
        let mut rotated = voxel_row(5, 1.0, 0.0, 0.1);
        let welded = boundary_edge_counts(&rotated, 5);
        assert_eq!(welded, alloc::vec![3, 2, 2, 2, 3]);

        // Dropping the tolerance can only lose cuts, never find new ones, and this rotation puts
        // at least one seam beyond an exact comparison.
        rotated.set_flags(CompoundFlags::FIX_INTERNAL_EDGES, Some(0.0));
        let exact = boundary_edge_counts(&rotated, 5);
        assert!(
            exact.iter().zip(welded.iter()).all(|(a, b)| a >= b),
            "exact found a cut the tolerance missed: {exact:?}"
        );
        assert!(
            exact.iter().sum::<usize>() > welded.iter().sum::<usize>(),
            "exact matched every seam: {exact:?}"
        );

        // A grid whose corners land on exact coordinates never needed the tolerance.
        let mut exact = voxel_row(5, 1.0, 0.0, 0.0);
        exact.set_flags(CompoundFlags::FIX_INTERNAL_EDGES, Some(0.0));
        assert_eq!(boundary_edge_counts(&exact, 5), alloc::vec![3, 2, 2, 2, 3]);
    }

    /// Every interior seam of a full grid is found, and only the perimeter survives: an `n` by `n`
    /// grid of cells encloses exactly `4 * n` boundary edges. Four cells meet at every interior
    /// corner, so this also exercises welding more than two coincident corners at once.
    #[test]
    fn a_voxel_grid_keeps_only_its_perimeter() {
        for &(n, size, origin, angle) in &[
            (3u32, 1.0 as Real, 0.0 as Real, 0.1 as Real),
            (3, 0.1, 0.0, 0.0),
            (4, 0.7, 5.0, 0.7854),
        ] {
            let half = size / 2.0;
            let rotation = Rotation::from_angle(angle);
            let parts = (0..n)
                .flat_map(|i| {
                    (0..n).map(move |j| {
                        let center =
                            Vector::new(origin + i as Real * size, origin + j as Real * size);
                        (
                            Pose::from_parts(rotation * center, rotation),
                            SharedShape::new(Cuboid::new(Vector::new(half, half))),
                        )
                    })
                })
                .collect();

            let mut compound = Compound::new(parts);
            compound.set_flags(CompoundFlags::FIX_INTERNAL_EDGES, None);

            let counts = boundary_edge_counts(&compound, n * n);
            assert_eq!(
                counts.iter().sum::<usize>(),
                (4 * n) as usize,
                "size {size} origin {origin} angle {angle}: {counts:?}"
            );
        }
    }

    /// The two signed zeros are the same point and hash apart, so a corner one part spells `0.0`
    /// and another `-0.0` has to weld even where no tolerance is allowed at all.
    #[test]
    fn signed_zero_corners_weld_at_zero_tolerance() {
        assert_ne!(
            (0.0 as Real).to_bits(),
            (-0.0 as Real).to_bits(),
            "the two zeros are meant to differ bitwise"
        );

        // Two triangles sharing the edge from the origin to `(1, 0)`, the second spelling its end
        // of that edge with negative zeros.
        let mut compound = Compound::new(vec![
            (
                Pose::IDENTITY,
                SharedShape::new(Triangle::new(
                    Vector::new(0.0, 0.0),
                    Vector::new(1.0, 0.0),
                    Vector::new(0.5, 1.0),
                )),
            ),
            (
                Pose::IDENTITY,
                SharedShape::new(Triangle::new(
                    Vector::new(-0.0, -0.0),
                    Vector::new(0.5, -1.0),
                    Vector::new(1.0, -0.0),
                )),
            ),
        ]);
        compound.set_flags(CompoundFlags::FIX_INTERNAL_EDGES, Some(0.0));

        assert_eq!(boundary_edge_counts(&compound, 2), alloc::vec![2, 2]);
    }

    /// The tolerance follows each corner's own distance from the origin rather than the whole
    /// compound's, so a part placed far away does not coarsen detail near it.
    #[test]
    fn a_distant_part_does_not_coarsen_detail_near_the_origin() {
        let small = 1.0e-3;
        let box_at = |x: Real, half: Real| {
            (
                Pose::from_parts(Vector::new(x, 0.0), Default::default()),
                SharedShape::new(Cuboid::new(Vector::new(half, half))),
            )
        };

        // Two small boxes with a gap between them, and one part 1e5 away. Scaled to the distant
        // part the gap would vanish; scaled to the small boxes themselves it is enormous.
        let mut compound = Compound::new(vec![
            box_at(0.0, small),
            box_at(4.0 * small, small),
            box_at(1.0e5, 1.0),
        ]);
        compound.set_flags(CompoundFlags::FIX_INTERNAL_EDGES, None);

        for part in 0..2 {
            assert_eq!(
                compound
                    .part_normal_constraints(part)
                    .unwrap()
                    .boundary_edges
                    .len(),
                4,
                "part {part} lost an edge to the distant part's scale"
            );
        }
    }

    #[test]
    fn contacts_with_the_ramp_itself_are_left_alone() {
        let plain = ramp_foot();
        let mut fixed = ramp_foot();
        fixed.set_flags(CompoundFlags::FIX_INTERNAL_EDGES, None);

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
    use crate::math::{Pose, Real, Rotation, Vector};
    use crate::shape::{Cuboid, SharedShape};
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
        compound.set_flags(CompoundFlags::FIX_INTERNAL_EDGES, None);

        // Each prism has 6 faces and loses the shared quad -- which only holds if `convex_hull`
        // produced the cut with the same corners on both sides, the premise everything rests on.
        for part in 0..2 {
            let cones = compound.part_normal_constraints(part).expect("polyhedron");
            assert_eq!(cones.boundary_faces.len(), 5, "part {part}");
        }
    }

    /// The 2D cut L, extruded. The cut is a quad, and it reaches the union's surface along two
    /// vertical edges: the convex one at the origin, and the reflex one above `(1, 1)`.
    fn cut_ell() -> Compound {
        let prism = |outline: &[[Real; 2]]| {
            let corners: Vec<Vector> = outline
                .iter()
                .flat_map(|p| [Vector::new(p[0], p[1], -1.0), Vector::new(p[0], p[1], 1.0)])
                .collect();
            (
                Pose::IDENTITY,
                SharedShape::convex_hull(&corners).expect("convex prism"),
            )
        };

        let mut compound = Compound::new(vec![
            prism(&[[0.0, 0.0], [2.0, 0.0], [2.0, 1.0], [1.0, 1.0]]),
            prism(&[[0.0, 0.0], [1.0, 1.0], [1.0, 2.0], [0.0, 2.0]]),
        ]);
        compound.set_flags(CompoundFlags::FIX_INTERNAL_EDGES, None);
        compound
    }

    /// The pseudo-normal the cone of `part`'s face pointing along `face` carries at the vertical
    /// edge standing over `corner`.
    fn pseudo_normal_at_vertical_edge(
        compound: &Compound,
        part: u32,
        face: Vector,
        corner: [Real; 2],
    ) -> Vector {
        let (bottom, top) = (
            Vector::new(corner[0], corner[1], -1.0),
            Vector::new(corner[0], corner[1], 1.0),
        );

        let faces = compound.part_faces();
        let (_, corners) = faces[part as usize]
            .as_ref()
            .expect("polyhedron")
            .iter()
            .find(|(normal, _)| (*normal - face).length() < 1.0e-5)
            .unwrap_or_else(|| panic!("no face pointing along {face:?}"));

        // `edge_pseudo_normals[i]` describes the edge from `corners[i]` to its successor.
        let edge = (0..corners.len())
            .find(|i| {
                let ends = [corners[*i], corners[(*i + 1) % corners.len()]];
                ends.iter().any(|end| (*end - bottom).length() < 1.0e-5)
                    && ends.iter().any(|end| (*end - top).length() < 1.0e-5)
            })
            .unwrap_or_else(|| panic!("no vertical edge over {corner:?}"));

        compound
            .part_normal_constraints(part)
            .expect("polyhedron")
            .boundary_faces
            .iter()
            .find(|cone| (cone.face - face).length() < 1.0e-5)
            .expect("cone")
            .edge_pseudo_normals[edge]
    }

    /// An edge a cut splits is still an edge of the union, so the two faces meeting along it open
    /// towards each other even though they belong to different parts.
    #[test]
    fn a_convex_edge_split_by_a_cut_opens_across_the_parts() {
        let compound = cut_ell();

        // Each prism keeps five faces and loses the cut quad it shares.
        for part in 0..2 {
            assert_eq!(
                compound
                    .part_normal_constraints(part)
                    .unwrap()
                    .boundary_faces
                    .len(),
                5,
                "part {part}"
            );
        }

        // The union turns 90 degrees along the vertical edge over the origin, between the front
        // face of one prism and the left face of the other.
        let halfway = Vector::new(-1.0, -1.0, 0.0).normalize();
        for (part, face) in [(0, -Vector::Y), (1, -Vector::X)] {
            let pseudo_normal = pseudo_normal_at_vertical_edge(&compound, part, face, [0.0, 0.0]);
            assert!(
                (pseudo_normal - halfway).length() < 1.0e-5,
                "part {part} stopped at {pseudo_normal:?}"
            );
        }
    }

    /// A concave edge has no outward range for its faces to share, so each keeps its own normal --
    /// the same thing that happens along an edge whose other face was cut away.
    #[test]
    fn a_concave_edge_split_by_a_cut_stays_closed() {
        let compound = cut_ell();

        // The reflex edge over `(1, 1)`, where the step rises out of the lower arm of the L.
        for (part, face) in [(0, Vector::Y), (1, Vector::X)] {
            let pseudo_normal = pseudo_normal_at_vertical_edge(&compound, part, face, [1.0, 1.0]);
            assert!(
                (pseudo_normal - face).length() < 1.0e-5,
                "part {part} opened to {pseudo_normal:?}"
            );
        }
    }

    /// A row of cuboids on a shared grid, the way a voxel body reaches a compound, optionally
    /// rotated as a whole.
    fn voxel_row(n: u32, size: Real, origin: Real, angle: Real) -> Compound {
        let half = size / 2.0;
        let rotation = Rotation::from_rotation_z(angle);
        let parts = (0..n)
            .map(|i| {
                let center = Vector::new(origin + i as Real * size, 0.0, 0.0);
                (
                    Pose::from_parts(rotation * center, rotation),
                    SharedShape::new(Cuboid::new(Vector::new(half, half, half))),
                )
            })
            .collect();

        let mut compound = Compound::new(parts);
        compound.set_flags(CompoundFlags::FIX_INTERNAL_EDGES, None);
        compound
    }

    /// How many boundary faces each part of a row kept. A welded row leaves the two ends with
    /// five and every box between them with four.
    fn boundary_face_counts(compound: &Compound, n: u32) -> Vec<usize> {
        (0..n)
            .map(|i| {
                compound
                    .part_normal_constraints(i)
                    .map_or(99, |cones| cones.boundary_faces.len())
            })
            .collect()
    }

    /// A cuboid is a polyhedron the type system does not call one, so a voxel body used to get no
    /// cones at all; and its parts compute their own corners, which need not round alike.
    #[test]
    fn a_row_of_cuboids_shares_its_faces() {
        for &(size, origin, angle) in &[
            (1.0 as Real, 0.0 as Real, 0.0 as Real),
            (0.1, 0.0, 0.0),
            (0.7, 0.0, 0.0),
            (1.0, 0.05, 0.0),
            (0.1, 1000.0, 0.0),
            (1.0, 0.0, 0.1),
            (1.0, 1000.0, 0.7854),
        ] {
            assert_eq!(
                boundary_face_counts(&voxel_row(5, size, origin, angle), 5),
                vec![5, 4, 4, 4, 5],
                "size {size} origin {origin} angle {angle}"
            );
        }
    }

    /// A zero tolerance is the exact comparison, which only recognizes a cut whose corners round
    /// to the very same bits on both sides.
    #[test]
    fn a_zero_tolerance_welds_only_bit_identical_corners() {
        let mut rotated = voxel_row(5, 1.0, 0.0, 0.1);
        let welded = boundary_face_counts(&rotated, 5);
        assert_eq!(welded, vec![5, 4, 4, 4, 5]);

        // Dropping the tolerance can only lose cuts, never find new ones, and this rotation puts
        // at least one seam beyond an exact comparison.
        rotated.set_flags(CompoundFlags::FIX_INTERNAL_EDGES, Some(0.0));
        let exact = boundary_face_counts(&rotated, 5);
        assert!(
            exact.iter().zip(welded.iter()).all(|(a, b)| a >= b),
            "exact found a cut the tolerance missed: {exact:?}"
        );
        assert!(
            exact.iter().sum::<usize>() > welded.iter().sum::<usize>(),
            "exact matched every seam: {exact:?}"
        );

        // A grid whose corners land on exact coordinates never needed the tolerance.
        let mut exact = voxel_row(5, 1.0, 0.0, 0.0);
        exact.set_flags(CompoundFlags::FIX_INTERNAL_EDGES, Some(0.0));
        assert_eq!(boundary_face_counts(&exact, 5), vec![5, 4, 4, 4, 5]);
    }

    /// A body pressed into the wedge tip reports a normal along its own face, not a surface of
    /// the union; the cones pull it onto the one face that is.
    #[test]
    fn a_tip_normal_is_clamped_onto_the_ramp_face() {
        use crate::query::details::NormalConstraints;

        let mut compound = ramp_foot();
        compound.set_flags(CompoundFlags::FIX_INTERNAL_EDGES, None);
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
