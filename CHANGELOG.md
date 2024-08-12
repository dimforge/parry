# Change Log

## Unreleased

### Modified

- Improve convergence of epa algorithm in degenerate configurations.

## v0.17.0

### Added

- Add `Triangle::robust_scaled_normal` and `Triangle::robust_normal` as a more robust way to compute the triangles
  normal for thin triangles that generally cause numerical instabilities.
- Add `Triangle::angle_closest_to_90` to find the triangle’s vertex with an angle closest to 90 degree.
- Add the `wavefront` feature that enables `TriMesh::to_obj_file` for exporting a mesh as an obj file.
- Add `Shape::scale_dyn` for scaling a shape as a trait-object.

### Modified

- `TypedShape::Custom(u32)` is now `TypedShape::Custom(&dyn Shape)`.
- `AabbSetsInterferencesCollector::tolerence` is now spelled correctly as `tolerance`.
- `Real` is now exposed through a `use` statement,
  so that an indirection is removed in documentation:
  previous occurrences of `Real` now show `f32` or `f64`.
- Significantly improved the general stability of mesh/mesh intersection calculation.
- Rename `Shape::clone_box` to `Shape::clone_dyn` (the `clone_box` method still exists but has been
  deprecated).
- Make `try_convex_hull` return an error instead of panicking if less than 3 input points are given.
- Make `Triangle::normal` and `Triangle::scaled_normal` only available in 3D instead of panicking in 2D.

## v0.16.1

### Fix

- Fix occasional crash in mesh/mesh intersection if some of the vertex coordinates are very small.

## v0.16.0

### Fix

- Fix edge case where some of the principal angular inertia are clamped to zero
  for decimeter-sized objects.
- Have ball-ball shape casting take into account the `stop_on_penetration` flags.
- Don’t panic in EPA for a corner case that needs some additional debugging. Show a debug log instead.

### Added

- Implement concave polygons intersections: `polygons_intersection_points`, `polygon_intersection`.

### Modified

- Update `bitflags` to version ^2.3
- Update `nalgebra` to 0.33.
- Update `indexmap` to 2.

## v0.15.1

### Fix

- Fix a regression in ball vs. convex shape contact manifold calculation.

## v0.15.0

### Added

- Add `ShapeCastOptions` that includes two new options for (linear) shape-casting.
  `ShapeCastOptions::target_distance` which will return a hit as soon as the moving
  shapes are closer than this distance; and `compute_impact_geometry_on_penetration`
  which forces the calculation of proper witness points and normals even if the shapes
  are initially intersecting (`time_of_impact == 0.0`).

### Modified

This version modifies many names related to shape-casting:

- Renamed `TOI` to `ShapeCastHit`.
- Renamed `TOIStatus` to `ShapeCastStatus`.
- Rename `RayIntersection::toi` to `RayIntersection::time_of_impact`.
- More generally, all occurrences of the word `toi` have been replaced by `time_of_impact`
  for better clarity.
- Rename `query::time_of_impact` to `query::cast_shapes`. More generally, all the
  functions prefixed with `time_of_impact_` (e.g. `time_of_impact_ball_ball`) are
  now prefixed with `cast_shapes_` (e.g. `cast_shapes_ball_ball`).
- Rename `QueryDispatcher::time_of_impact` to `QueryDispatcher::cast_shapes`.
- The (linear) shape-casting functions like `query::cast_shapes` (previously named
  `query::time_of_impact) now take a `ShapeCastOptions` instead of the `max_toi` and
  `stop_at_penetration` arguments.
- Rename `query::nonlinear_time_of_impact` to `query::cast_shapes_nonlinear`.
- Rename `QueryDispatcher::nonlinear_time_of_impact` to `QueryDispatcher::cast_shapes_nonlinear`.
- Rename `NonlinearTOIMode` to `NonlinearShapeCastMode`, and `NonlinearTOIMode::DirectionalTOI` to
  `NonlinearShapeCastMode::Directional`.
- Rename `TimeOfImpactStatus::Penetrating` to `ShapeCastStatus::PenetratingOrWithinTargetDist`.

## v0.14.0

### Modified

- Remove CUDA support to break free from the toolchain restriction required by cust.
- Rework internal edges resolution using normal cones. This implies the modification of the
  `SimdCompositeShape::map_part_at`, `TypedSimdCompositeShape::map_typed_part`, and
  `TypedSimdCompositeShape::map_untyped_part` trait functions so that the closure argument takes
  an extra argument for the (optional) normal constraints. This argument can be safely ignored
  by user code unless applying the normal collection is relevant to your use-case.
- Contact manifolds will now retain all contacts (including the ones further than the specified `prediction`
  distance) whenever any contact is actually closer than this `prediction` distance.
- Typo fix: renamed `TopologyError::BadAdjascentTrianglesOrientation` to `BadAdjacentTrianglesOrientation`.

### Fixed

- Fix contacts between convex shapes being occasionally ignored due to some rounding errors.
- Remove crash when entering unreachable code in non-linear TOI calculation.
- Fix accuracy issue in triangle-mesh center-of-mass calculation when the mesh isn’t manifold.

### Added

- Add `SdpMatrix2::inverse_and_get_determinant_unchecked`. This is useful for computing the
  inverse in a AoSoA SIMD setting.
- Add `Aabb::intersects_moving_aabb` to perform a swept test between two moving aabbs.

## v0.13.8

### Added

- Add `Qbvh::traverse_depth_first_with_context`,  `Qbvh::traverse_depth_first_node_with_stack_and_context`, and the
  related `SimdVisitorWithContext` trait to allow parent nodes to pass a custom context to its children during
  recursion.

## v0.13.7

### Modified

- The `point_in_poly2d` now handles arbitrary (convex and non-convex) polygons. The previous implementation
  that only supported convex polygons has been renamed `point_in_convex_poly2d`.

### Fixed

- Fix a crash in `Qbvh::refit` that results from the QBVH tree becoming increasingly imbalanced.

### Added

- Add `Aabb::scaled_wrt_center` to scale an AABB while keeping its center unchanged.

## v0.13.6

### Fixed

- Fix ball-convex manifolds missing contacts in some corner cases.
- Fix panic in `TriMesh::intersection_with_plane`.

## v0.13.5

### Fixed

- When using `rkyv`, fix `CheckBytes` implementation for types archived as themselves.
- Fix occasional crash in the `QBVH` incremental update.

## v0.13.4

### Fixed

- Fix `Polyline::flat_indices` that returned an incorrectly sized slice.
- Fix serialization of `SimdAabb` into map-like formats like JSON, YAML, RON.

### Added

- Add validation when using `rkyv` safe API whenever applicable to `parry` types.

## v0.13.3 (08 March 2023)

### Modified

- Improved performance of intersection checks involving composite shapes (compound shapes, trimeshes, polylines, etc.)

## v0.13.2 (08 March 2023)

This version was yanked. See the release notes for 0.13.3 instead.

## v0.13.1 (26 Feb. 2023)

### Fixed

- Add workaround to address jitter issue due to incorrectly empty contact manifolds generated sometimes for convex  
  polyhedron.

## v0.13.0 (15 Jan. 2023)

### Modified

- About `rkyv` support: most POD structs (`Aabb`, `Ball`, `Cuboid`, etc.) are now archived as themselves instead of
  being archived as different types (for example `Aabb` is archived as `Aabb` itself instead of `ArchivedAabb`).

### Added

- In 3D, add `transformation::try_convex_hull` for a convex hull calculation that will return an error instead of
  panicking on unsupported inputs.

### Fixed

- Fixed duplicate faces in the connected components returned by `TriMesh::connected_components`.

## v0.12.1 (09 Jan. 2023)

### Added

- Add `TriMesh::canonical_intersection_with_plane` for intersecting with planes aligned with one of the coordinate axes.
- Add `TriMesh::intersection_with_plane` for intersecting with arbitrary planes.
- Add `TriMesh::intersection_with_local_plane` for intersecting with arbitrary planes in the same space as the mesh
- Add `IntersectResult` as the output type for the above functions.
- Add `Polyline::extract_connected_components` which splits a compound polyline into its connected components.
- Add implementations of `bytemuck::Pod` and `bytemuck::Zeroable` for all the simple shapes that allow it
  (`Cuboid`, `Ball`, `Cone`, etc.), and for bounding volumes (`BoundingSphere` and `Aabb`).

## v0.12.0 (11 Dec. 2022)

### Modified

- `Qbvh::leaf_data` now requires `&self` instead of `&mut self`.
- Replace the `Qbvh::leaf` boolean by a bitflags.

### Added

- Add `Qbvh::remove`, `Qbvh::pre_update_or_insert`, `Qbvh::refit`, `Qbvh::rebalance` to allow modifying a `Qbvh`
  without having to rebuild it completely.
- Add `QbvhNode::is_leaf` to get if a node is a leaf or not.
- Add `SharedShape::trimesh_with_flags` for building a trimesh with specific pre-processing flags.

### Fixed

- Fix `Triangle::contains_point`.

## v0.11.1 (30 Oct. 2022)

### Added

- Add `SharedShape::trimesh_with_flags` for constructing a triangle mesh with flags
  specified by the user.

## v0.11.0 (30 Oct. 2022)

### Modified

- Rename `AABB` to `Aabb` to comply with Rust’s style guide.
- Rename `QBVH` to `Qbvh` to comply with Rust’s style guide.

### Added

- Add `ConvexPolygon::offsetted` to dilate a polygon.
- Add `CudaTriMesh` and `CudaTriMeshPtr` for triangle-meshes usable with CUDA.
- Add a no-std implementation of point-projection on a triangle mesh.

### Fixed

- Fix ghost collisions on internal edges on flat 3D meshed and flat 3D heightfields.
- Fix pseudo-normals calculation that could generate invalid normals for triangles with
  some small vertex angles.
- Fix `Aabb::bounding_sphere` which returned a bounding sphere that was too big.

## v0.10.0 (02 Oct. 2022)

### Modified

- Add to `query::time_of_impact` a boolean argument `stop_at_penetration`. If set to `false`
  the linear shape-cast won’t immediately stop if the shape is penetrating another shape at its
  starting point **and** its trajectory is such that it’s existing that penetration configuration.

### Added

- Add 2D `Heightfield::to_polyline` to get the explicit vertices/indices of a 2D heightfield
  seen as a polyline.
- Add the support for linear shape-cast (`query::time_of_impact`) for heightfields.
- Make the convex polyhedron scaling more forgiving regarding normals to avoid frequent unjustified panics.
- Fix panic happening when building a convex polyhedron with empty inputs.
- Add the support of Heightfields on CUDA kernels written in Rust using the `cust` crate.
- Add the `rkyv-serialize` feature that enables the implementation of `rkyv` serialization/deserialization
  for most shapes.
- Add the `parallel` feature that enables methods for the parallel traversal of Qbvh
  trees: `Qbvh::traverse_bvtt_parallel`,
  `Qbvh::traverse_bvtt_node_parallel`, `Qbvh::traverse_depth_first_parallel`, `Qbvh::traverse_depth_first_node_parallel`.

### Fixed

- Fix the application of non-uniform scaling to balls.

## v0.9.0 (30 Apr. 2022)

### Modified

- Remove `&self` argument from `Compound::decompose_trimesh`.
- Switch to `cust` 0.3 (for partial CUDA support).
- Rename `RoundShape::base_shape` to `RoundShape::inner_shape`.

### Added

- Allow custom balancing strategies for the Qbvh construction. Some strategies are allowed to generate
  new leaves during the splitting process.
- Allow using point projection on heightfields from a CUDA kernel.
- Add the simultaneous traversal of two Qbvhs.
- Add computation of `MassProperties` for a `TriMesh`.
- Add `.to_outline` methods to compute the outline of a 3D shape (useful for debug-rendering).
- Add method to apply a scaling factor to some shapes. Shapes not supporting non-uniform scaling (like balls)
  will return a convex approximation of the scaled result.
- Add methods to split (into up to two parts) a Cuboid, Segment, or TriMesh with an axis-aligned plane.
- Add the computation of the intersection mesh between two TriMesh, or between a Cuboid and a TriMesh.

## v0.8.0 (2 Jan. 2022)

### Modified

- Until now, the orientation of the polygon computed by 2D convex hull
  computation `parry2d::transformation::convex_hull` and
  `parry2d::transformation::convex_hull_idx` wasn't specified (and was generally in clockwise order). Now, this
  orientation is explicitly specified in the documentation and is set to counter-clockwise order (which is coherent
  with orientation expected by, e.g., the `ConvexPolygon` type).

### Added

- Add `parry::utils::obb` which computes a (possibly sub-optimal) OBB for a set of points.
- Add `Polyline::project_local_point_assuming_solid_interior_ccw` which projects a point on the polyline contour, and
  is able to detect if that points is located inside of the polyline, assuming that the polyline is closed and oriented
  counter-clock-wise.
- Add (3D only) `TriMesh::compute_pseudo_normals`. If this is called, and if the trimesh respects some constraints (
  oriented
  with outward normals, manifold almost everywhere, etc.) any subsequent point queries with the `solid` argument
  set to `true` will properly set the `PointProjection::is_inside` to true when the point lies in the interior of
  the trimesh.
- Added the implementation of the ear-clipping and Hertel-Mehlhorn algorithm that can be used for 2D triangulation and
  convex decomposition.
- Add the ability to use a subset of Parry’s features in a `no-std` context. If the new `cuda` cargo feature of Parry is
  enabled, all features compatible with `no-std` can be used inside of a CUDA kernel written in Rust thanks to
  the [rust-cuda](https://github.com/Rust-GPU/Rust-CUDA) ecosystem.

### Fixed

- Fix the orientation of the polygons generated by the 2D convex polygon decomposition (they are now always oriented
  counter-clockwise as expected by the `ConvexPolygon` type).
- Fix the intersection test between a 2D ball and a 2D compound shape.

## v0.7.1

### Added

- Add the method `Aabb::volume` to compute the volume of an Aabb.

## v0.7.0

### Modified

- Update the codebase to use `nalgebra v0.29`.

### Fixed

- Fix a bug where the normal returned by ray-casting on polylines would not be normalized.

## v0.6.0

### Added

- Implement `Debug, Clone, PartialEq` for `VHACDParameters`.
- Add a method to reverse the order of a polyline.
- Add a method to remove duplicate vertices form a `TriMesh` (and adjusting the index buffer accordingly).
- Add a method to iterate through all the lean data stored by a Qbvh.
- Implement the Interval Newton Method for computing all the roots of a non-linear scalar function.
- Implement the intersection test between a spiral and an Aabb.

### Modified

- Rename all occurrences of `quadtree` to `qbvh`. Using the term `quadtree` was not representative of the actual
  acceleration structure being used (which is a BVH).

## v0.5.1

### Fixed

- Fix a bug where `query::contact` would return `None` for a all intersecting a cuboid
  in such a way that one of the cuboid's vertices coincides exactly with the ball's
  center.

## v0.5.0

### Modified

- Updated all dependencies to their latest version.

### Fixed

- Fix ray-casting against solid triangles.
- Fix NaN when adding mass-properties with a zero mass.

## v0.4.2

### Added

- `ShapeType` now implements `PartialEq, Eq, Hash`.

### Fixed

- The order of vertices output by `Cuboid::to_polyline` has been modified
  to actually represent the cuboid's boundary (instead of passing through
  its diagonal).

## v0.4.1

### Added

- `SharedShape` now implements `AsRef<dyn Shape>`.
- Add the optional method `Shape::compute_swept_aabb` to the `Shape` trait.

### Modified

- Renamed `SimdQuadTree` to `Qbvh` (Quaternary Bounding Volume Hierarchy). The
  incorrect name `SimdQuadTree` is now deprecated.
- `Qbvh::clear_and_rebuild` is now slightly more general.

## v0.4.0

### Modified

- Switch to `nalgebra` 0.26.

## v0.3.0

### Added

- Add a special case for the triangle/cuboid edge-edge case in the SAT implementation.
- Add contact manifold computation between a convex shape and a HalfSpace.
- Add a special case (instead of GJK) for computing the distance between a ball and a convex shape.
- Add a `Shape::clone_box` method to clone a shape trait-object.
- Add a `Shape::make_mut` method to take a mutable reference to the shape trait-object.
- Add methods like `Shape::as_shape_mut` to downcast a `&mut Shape` to a concrete mutable shape.
- Add `MassProperties::set_mass` as a simple way to modify the mass while automatically adjusting
  the angular inertia.
- Make the fields of `BoundingSphere` public.
- Add `Aabb::EDGES_VERTEX_IDS` and `Aabb::FACES_VERTEX_IDS` which are index tables describing
  the vertices of a given edge, or of a given face, of the Aabb.

### Modified

- Remove the `target_dist` argument from `query::time_of_impact`.
- The `RigidBodyMotion` and all its implementors have been removed. Use the `NonlinearRigidMotion`
  struct instead.
- The `query::nonlinear_time_of_impact` function arguments changed completely to accommodate for
  the new `NonlinearRigidMotion` struct. See the doc of `query::nonlinear_time_of_impact` for details.

### Fixed

- Fix the separation computation between a cuboid and a support-map shape in the SAT implementation.
- Fix a bug causing invalid results with `query::time_of_impact` whenever the first shape was rotated.

## v0.2.0

### Added

- Add the `Shape::as_typed_shape` method to convert the shape trait-object
  to an enum.
- Add the `TypedShape` enum with variants containing references to concrete shapes
  (except for user-defined custom shapes).
- Add dedicated algorithms for the projection of a point on a cone or a cylinder.
- Improve the overall shape serialization performance.
- Add `Aabb::vertices` to get an array containing its vertices.
- Add `Aabb::split_at_center` to split an Aabb into 4 (in 2D) or 8 (in 3D) parts.
- Add `Aabb::to_trimesh` to compute the mesh representation of an Aabb.

### Removed

- Remove the `Shape::as_serialize` method.

### Fixed

- Fix a bug causing making some ball/convex shape contact manifold computation
  fail when they are penetrating deeply.
