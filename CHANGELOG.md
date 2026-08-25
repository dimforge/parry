## Unreleased

### Added

- `CompoundFlags::FIX_INTERNAL_EDGES` (2D) makes a `Compound` treat the edges its parts share as
  interior to the union, so a body sliding across the cut between two parts of a convex
  decomposition no longer catches on it. `Compound::PartNormalConstraints` is now
  `CompoundPseudoNormals`, matching what `TriMesh` and `Polyline` already provide.
- Ray casts and point projections against a `Compound` now identify the part that was hit:
  the returned feature is `FeatureId::Face(part_index)`, the way `TriMesh` reports the
  triangle. Both previously discarded the part index.

## 0.30.2

### Added

- `Debug` for `TriMesh` now prints a summary instead of just the type name
  ([#177](https://github.com/dimforge/parry/issues/177)).
- `query::contact` now supports heightfields, instead of returning `Err(Unsupported)`
  ([#274](https://github.com/dimforge/parry/issues/274)).

### Fixed

- Fix GJK reporting overlapping axisymmetric shapes as disjoint when the simplex stalls on a
  degenerate support direction ([#396](https://github.com/dimforge/parry/issues/396)).
- Fix wrong shape-cast normals and distances for stalled GJK simplices, now classified as touching
  or penetrating from the solver's own distance bound ([#193](https://github.com/dimforge/parry/issues/193)).
- Fix shape-casts missing grazing hits at large coordinates, by scaling the GJK ray-cast tolerance
  with the support magnitude ([#180](https://github.com/dimforge/parry/issues/180)).
- Fix shape-cast time of impact being short by an amount scaling with the shape's extent
  ([#429](https://github.com/dimforge/parry/issues/429)).
- Fix `ShapeCastStatus::PenetratingOrWithinTargetDist` being reported for merely touching support-map
  casts ([#106](https://github.com/dimforge/parry/issues/106)).
- Fix `contact` between exactly-touching shapes returning witness points at the shapes' centers
  ([#315](https://github.com/dimforge/parry/issues/315)).
- Fix `Triangle::area` returning wildly wrong values for degenerate triangles
  ([#111](https://github.com/dimforge/parry/issues/111)).
- Fix point projection on degenerate triangles reporting the point as inside
  ([#76](https://github.com/dimforge/parry/issues/76)).
- Fix `segments_intersection2d` classifying exact segment endpoints as `OnEdge` instead of `OnVertex`
  ([#109](https://github.com/dimforge/parry/issues/109)).
- Fix composite-shape point projection panicking on non-finite points
  ([#395](https://github.com/dimforge/parry/issues/395)).
- Fix `HeightField::project_local_point` iterating on every element, making it O(rows * cols)
  ([rapier#332](https://github.com/dimforge/rapier/issues/332)).
- Fix single-point contact manifolds on cylinder and cone caps, by orienting the cap's polygonal
  approximation toward the contact ([rapier#810](https://github.com/dimforge/rapier/issues/810)).
- Fix missing speculative contacts between a `Voxels` shape and a shape within the prediction
  distance ([#404](https://github.com/dimforge/parry/issues/404)).
- Fix `Voxels` shape-cast `witness1` not being expressed in the voxels shape's local frame
  ([#373](https://github.com/dimforge/parry/issues/373)).
- Fix an out-of-bounds panic in the binned BVH builder on degenerate leaf AABBs
  ([rapier#961](https://github.com/dimforge/rapier/issues/961)).
- Fix out-of-bounds indices in `to_outline` for round shapes with a zero border radius
  ([rapier#969](https://github.com/dimforge/rapier/issues/969)).
- Fix the quadratic polygon removal in `hertel_mehlhorn`
  ([#408](https://github.com/dimforge/parry/issues/408)).
- Require `spade` 2.15, the minimum version providing `try_bulk_load_cdt`
  ([#428](https://github.com/dimforge/parry/issues/428)).

### Modified

- Document the frame conventions of the `ShapeCastHit` fields
  ([rapier#933](https://github.com/dimforge/rapier/issues/933)).
- Document that the `Voxels` iterators only yield non-empty voxels
  ([#382](https://github.com/dimforge/parry/issues/382)).
- Document that `BoundingVolume::merge`/`merged` don't guarantee strict containment of their inputs
  ([#260](https://github.com/dimforge/parry/issues/260)).

## 0.30.1

### Fixed

- Fix incorrect `closest_points_cuboid_cuboid` closest points (affecting `query::distance` too) in degenerate
  configurations where the separating axis is ambiguous, e.g., axis-aligned cuboids closest at a
  vertex or along parallel edges.

## 0.30.0

### Breaking changes

- The `simd-stable` and `simd-nightly` features were removed. 4-lane SIMD is now always enabled
  (it falls back to scalar code on targets without SIMD support). The new opt-in `simd8` feature
  widens SIMD to 8 lanes for `f32` builds; it requires an AVX-enabled target to actually emit
  256-bit instructions, and is incompatible with `enhanced-determinism`.
- `ContactManifold::subshape_pos1`/`subshape_pos2` are no longer public fields. They are replaced by
  a single boxed `subshape_poses: Option<Box<SubshapePoses>>` field (to shrink `ContactManifold`),
  accessed through the new `subshape_pos1()`/`subshape_pos2()` getters and
  `set_subshape_pos1`/`set_subshape_pos2` setters.

### Added

- `query::sweep_toi`: sweep-based time-of-impact queries. A timestep is modeled as a `Sweep` between
  two endpoint poses (linear translation + rotation nlerp), and `sweep_time_of_impact` (plus
  `sweep_time_of_impact_composite` for composite shapes) computes the earliest time the swept shapes
  reach a slop-based target separation, using conservative advancement. Also exports `ToiProxy`,
  `SimplexCache`, `SweepToiOutput`, `SweepToiStatus`, and `SweepCompositeFastShape`.
- `Bvh` gains incremental and parallel update APIs: `refit_partial`, flag-preserving
  `refit_without_resolve` variants, `refit_parallel`, parallel BVTT traversal, and batched
  parallel leaf updates (`insert_or_update_batch_partially_parallel`,
  `reinsert_or_update_with_change_detection`, `reinsert_or_update_if_present`,
  `update_partially_if_present`, and the `BvhLeafUpdateStatus` enum).

### Modified

- Cuboid support-face feature ids are now computed with bit operations instead of lookup tables.

## 0.29.0

### Breaking changes

- The `<Polyline as TypedCompositeShape>::PartNormalConstraints` is now `SegmentPseudoNormals` instead of `()`.

### Fixed

- Fix `enhanced-determinism` not producing cross-platform deterministic results for parry3d/parry3d-f64.
  The feature now enables `glamx/scalar-math` to disable architecture-specific SIMD (NEON on arm64,
  SSE2 on x86_64) which caused floating-point non-associativity in Vec3/Vec4/Quat dot products.
- Oriented (one-sided) 2D polylines ([#254](https://github.com/dimforge/parry/issues/254)). Build a `Polyline` with
  `PolylineFlags::ORIENTED` via `Polyline::with_flags`, or toggle it with `Polyline::set_flags`, to clamp contact normals
  to the outward side. This removes the spurious sideways push a body gets at a convex corner of a thin double-sided
  polyline. A single 2D flag covers what `TriMesh` splits across `TriMeshFlags::ORIENTED` (assume an outward orientation
  and compute pseudo-normals) and `TriMeshFlags::FIX_INTERNAL_EDGES` (use them to clamp contact normals). Also adds
  `Polyline::flags`, `Polyline::segment_normal_constraints`, and the `PolylineFlags` / `SegmentPseudoNormals` types.
  ([#425](https://github.com/dimforge/parry/pull/425))

## 0.28.0

### Breaking changes

- `CompositeShapeRef::project_local_point` and `project_local_point_and_get_feature` now return
  `Option<...>` instead of unwrapping internally, and take an additional `max_dist` parameter that
  bounds the search distance.

### Modified

- Bump `glamx` to 0.3 (built on `glam` 0.33), `simba` to 0.10, `rstar` to 0.13, `hashbrown` to 0.17,
  `rand` to 0.10, and `kiss3d` (visual examples) to 0.42.
- Many shape, AABB, and voxel internals were reworked to avoid bracket-indexing of `glam` vectors,
  which is not supported on SPIR-V targets. A new `VectorExt::vget`/`vset` API and a `for_each_dim!`
  macro are now used in place of `v[i]`. `parry3d`'s `alloc` feature also now gates `smallvec`,
  `downcast-rs`, `rstar`, and `glamx/approx`, so the no-alloc build can target SPIR-V.

### Fixed

- Fix `Bvh::scale`/`BvhNode::scale` producing an inverted AABB on negative scale components, which
  silently broke queries on shapes scaled with reflections (e.g. `TriMesh::scaled`). Scale components
  must now be non-zero (previously: positive) ([#406](https://github.com/dimforge/parry/pull/406)).
- Fix `Voxels::scaled` not scaling the internal chunk BVH, causing queries on scaled voxels to miss
  geometry ([#419](https://github.com/dimforge/parry/pull/419)).
- Fix `Bvh::optimize_incremental` panicking after multiple `remove`s without an intervening `refit`
  when the root collapses to a partial root, leaving orphaned wide nodes
  ([#409](https://github.com/dimforge/parry/pull/409)).
- Correct the assertion message on `Aabb`/`BoundingSphere`'s `loosen`/`tighten` (and their `*_ed`
  variants) to "must be non-negative" ([#394](https://github.com/dimforge/parry/pull/394)).
- Fix multiple EPA failures on large coordinate magnitudes by scaling the face-rejection tolerance
  relative to the simplex vertex magnitudes ([#415](https://github.com/dimforge/parry/issues/415)).
- Fix CCD edge cases (in both ball-vs-ball and support-map-vs-support-map casts) where casts
  starting at (or very close to) the contact boundary with a near-tangent direction would return
  an unreliable normal. The fallback now uses the contact query to recover a robust closest-point
  normal, and the small-TOI threshold was relaxed from `1e-5` to `1e-4`.
- Fix `Aabb::cast_local_ray_and_get_normal` panicking with a subtraction overflow when casting a
  zero-direction ray starting inside the AABB ([#383](https://github.com/dimforge/parry/issues/383)).
- Fix infinite loop in `TriMesh::intersection_with_local_plane` on degenerate adjacency graphs
  that don't loop cleanly back to the starting index
  ([#398](https://github.com/dimforge/parry/issues/398)).
- Fix panic in `mesh_intersection` when constraint edges overlap (e.g. for co-planar triangles).
  Overlapping constraints are now skipped instead of crashing the CDT
  ([#389](https://github.com/dimforge/parry/issues/389)).
- Fix panic in 3D voxelization's internal convex-hull step by falling back to `try_convex_hull`
  and returning an empty hull on failure ([#347](https://github.com/dimforge/parry/issues/347)).
- Fix `WSign::copy_sign_to` on non-CUDA targets by using the native `copysign` (the bit-twiddling
  workaround is now scoped to `nvptx64`, where `cuda_std`'s `copysign` does not compile).

## 0.27.0

- Switched to `glamx 0.2` (`glam 0.32`).

## 0.26.1

### Added

- Added `Bvh::traverse_indexed` to traverse only a part of the subtree, and to retrieve the node’s indices ([#416](https://github.com/dimforge/parry/pull/416)).

### Fixed

- Fix `Bvh::optimize_incremental` being non-deterministic when restored from a serialized snapshot. The
  incremental optimization state is now stored on the `Bvh` itself (instead of the non-serializable
  `BvhWorkspace`) so it survives serialization/deserialization ([#403](https://github.com/dimforge/parry/pull/403)).

## 0.26.0

### Breaking changes

This release migrates parry from `nalgebra` to `glam` (via the `glamx` crate) for future compatibility with
`rust-gpu`. This is a major breaking change affecting almost all public APIs.

#### Type renames

- `Isometry` → `Pose` (using `glamx::Pose2`/`Pose3`)
- `Rotation` → `Rot2`/`Rot3` (using `glamx::Rot2`/`Rot3`)
- `IsometryOps` → `PoseOps`
- `IsometryOpt` → `PoseOpt`

#### Removed types

- `Point<Real>` - use `Vector` instead. Points and vectors are now unified.
- `UnitVector<Real>` - use `Vector` instead. Normalization is no longer encoded in the type.
- `Translation<Real>` - use `Vector` for translations.

#### Math type changes

- `Vector` is now `glam::Vec2`/`Vec3`/`DVec2`/`DVec3` depending on dimension and precision features
- `Matrix` is now `glam::Mat2`/`Mat3`/`DMat2`/`DMat3`
- The `math` module now re-exports glam types and provides dimension-agnostic aliases

#### API signature changes

- Many functions that previously took `&Point<Real>` or `&Vector<Real>` now take `Vector` by value
- Functions taking `&Isometry<Real>` now take `&Pose` or `Pose` by value
- `HalfSpace::new` now takes `Vector` instead of `Unit<Vector<Real>>`
- Shape constructors like `Segment::new`, `Triangle::new`, `Capsule::new` now take `Vector` instead of `Point`
- `Aabb::mins` and `Aabb::maxs` are now `Vector` instead of `Point<Real>`

#### Migration guide

If your codebase currently relies on `nalgebra`, note that `nalgebra` and `glamx` provide type conversion. Enable the
corresponding features:
- `nalgebra = { version = "0.34", features = [ "convert-glam030" ] }`
- `glamx = { version = "0.1", features = ["nalgebra"] }`
  then you can convert between `glam` and `nalgebra` types using `.into()`.

```rust
// Before (nalgebra)
use parry3d::na::{Point3, Vector3, Isometry3, Unit};
let point = Point3::new(1.0, 2.0, 3.0);
let vector = Vector3::new(1.0, 0.0, 0.0);
let normal = Unit::new_normalize(vector);
let pos = Isometry3::translation(1.0, 2.0, 3.0);

// After (glam)
use parry3d::math::{Vector, Pose};
let point = Vector::new(1.0, 2.0, 3.0);  // Points are now Vector
let vector = Vector::X;
let normal = vector.normalize();          // No Unit wrapper
let pos = Pose::translation(1.0, 2.0, 3.0);
```

Common patterns:
- `Point3::origin()` → `Vector::ZERO`
- `Vector3::x()` → `Vector::X`
- `point.coords` → just use the vector directly
- `Unit::new_normalize(v)` → `v.normalize()`
- `Isometry3::identity()` → `Pose::IDENTITY`
- `isometry.translation.vector` → `pose.translation`

### Modified

- Re-export `glamx` instead of `nalgebra` as the public linear algebra dependency
- Migrate visual examples from `macroquad` to `kiss3d`
- Renamed `rkyv-serialize` feature to just `rkyv`

## 0.25.3

- Significantly improve performances of `Voxels::combine_voxel_states`.

## 0.25.2

- Add the `encase` feature to enable implementation of `encase::ShaderType` for primitive shapes.
- Fix regression in shape-vs-composite-shape contact calculation.

## 0.25.1

- Fix point-projection of a ball on a voxels shape containing empty chunks (e.g. after voxels removal).
- General documentation improvements.

## 0.25.0

- The `Voxels` shape now uses a sparse storage internally.
- Remove `Voxels::try_set_voxel`. Use `::set_voxel` instead.
- The `Voxels::domain` function returns values instead of references now.
- Removed `Voxels::extents` and `Voxels::domain_center`. They can be obtained
  from `Voxels::local_aabb`.

## 0.24.0

- Fixed same-machine determinism of the re-exported hashmap when building without the `enhanced-determinism` feature.
- Replace all the square-root angular inertia tensors by the same tensor without taking it square root.
  As a result `MassProperties::world_inv_inertia_sqrt` has been renamed to `world_inv_inertia`,
  `MassProperties::inv_principal_inertia_sqrt` renamed to `inv_principal_inertia`.
- Modified the type aliases from the `simd` module when neither the `simd-stable` nor `simd-nightly` are enabled.
  Without these simd features all type aliases are single-lane non-simd types (`f32`, `f64`, `bool`) instead of 4-lanes auto-simd
  types.
- Moved the `SimdAabb` type behind the `simd-stable` or `simd-nightly` cargo feature.

## 0.23.0

- Update to nalgebra 0.34.

## 0.22.1

### Added

- Add `Bvh::insert_with_change_detection` which runs an insertion unless the inserted aabb is enclosed in the
  existing one. This differs from `insert_or_update_partially` because the former will update parent AABBs so
  they enclose their descendants.

## 0.22.0

### Fixed

- Fix bug in BVH tree node removal.
- Fix bug in BVH tree building from exactly two leaves.

## 0.22.0-beta.1

### Fixed

- Fix invalid BVH state that could be reached after a removal resulting in a partial root.

## 0.22.0-beta.0

### Added

- Add a new `Bvh` acceleration data-structure. It replaces `Qbvh` entirely. It supports:
  - Traversals (best-first, depth-first, BVTT, leaf iterators, and leaf pairs iterator).
  - It can be constructed either incrementally by inserting nodes, or from a set of leaves using either the
    binned building strategy or the PLOC (without parallelism) strategy.
  - Dynamic leaf insertion, update, removal.
  - Incremental tree rebalancing.

### Fixed

- Fix `clip_aabb_line` crashing when given incorrect inputs (zero length direction or NAN).
- Fix `Segment::intersects_ray` returning false-positive when the segment is zero-length. ([#31](https://github.com/dimforge/parry/issues/31)).
- Expose `utils::sort3` and `utils::sort2`.

### Modified

- The `local_point_cloud_aabb`, `point_cloud_aabb`, and `Aabb::from_points` now takes an iterator over point values instead
  of an iterator to point references. Variants taking point references still exist and are named `local_point_cloud_aabb_ref`,
  `point_cloud_aabb_ref` and `Aabb::from_points_ref`.
- Renamed `SimdCompositeShape` and `TypedSimdCompositeShape` to `CompositeShape` and TypedCompositeShape`.
- The `TypedCompositeShape` trait now derives from `CompositeShape`.
- Removed every `*Visitor` structures. Instead, either call `Bvh::traverse` (or `Bvh::search_best`, or `Bvh::leaves`, or
  `bvh::leaf_pairs`), or wrap your composite shape into `CompositeShapeRef` to access some generic implementation of
  various geometric queries for any composite shape.
- All composite shapes now rely on the new `Bvh` acceleration structure instead of `Qbvh`.
- The `Qbvh` has been removed. Use `Bvh` instead.

## 0.21.1

### Added

- Add `Voxels::combine_voxel_states` to merge the voxel state of two `Voxels` shapes, updating their internal voxels
  state as if both shapes were part of a single one. In particular, this will prevent any internal edge that would
  arise at the boundaries of both shapes if they were adjacent.
- Add `Voxels::propagate_voxel_change` to propagate a single-voxel modification from one `Voxels` shape to another,
  in order to update their internal neighborhood information as if both were part of the same `Voxels` shape. 

## 0.21.0

### Added

- Add `TriMesh::cast_ray_with_culling` and `TriMesh::cast_local_ray_with_culling` for casting rays on a triangle mesh
  but with the possibility to prevents hits on front-faces or back-faces.
- Add contact manifold calculation between two voxels shapes, or a voxels shape and compound shape.
- Add intersection check between voxels and other shapes.
- Add `MassProperties::from_voxels` to compute the mass and angular inertia tensor from a voxels shape.

### Modified

- Add new variants to `TypedWorkspaceData` for voxels-voxels and voxels-compound collision detection workspace
  data.
- The `Voxels` shape now only support cuboids as their leaf geometries (pseudo-balls were removed). 

## v0.20.2

### Fixed

- Fix infinite loop in `Voxels::set_voxel`.

## v0.20.1

### Added

- Rework the `Voxels` shape API to use better method names.
- Added implementations for linear and non-linear shape-casting involving `Voxels` shapes.

## v0.20.0 (yanked)

### Added

- Added the `Voxels` type: a dedicated shape for voxel models. This is currently experimental because some features are
  still missing (in particular: shape-casting, mass properties, and collision-detection against non-convex shapes).
- Added `SharedShape::voxels`, `SharedShape::voxels_from_points`, and `::voxelized_mesh` for creating a voxels shape
  from grid coordinates, points, or automatic voxelization of a triangle mesh.

## v0.19.0

### Added

- Derive `Copy` for `VHACDParameters`.
- Add `spade` default feature for algorithms using Delaunay triangulation from `spade`.
- Add `SharedShape::from_convex_polyline_unmodified` and `ConvexPolygon::from_convex_polyline_unmodified`
  to initialize a polyline from a set of points assumed to be convex, and without modifying this set even
  if some points are collinear.
- Add `TriMesh::pseudo_normals_if_oriented` that returns `Some` only if the mesh has the `TriMeshFlags::ORIENTED`
  flag enabled.

### Modified

- The `TriMeshFlags::FIX_INTERNAL_EDGES` flag no longer automatically enable the `TriMeshFlags::ORIENTED`
  flag (but the mesh pseudo-normals will still be computed).
- Improve `no_std` compatibility.
  - Everything is now compatible, except `mesh_intersections`, `split_trimesh`,
    convex hull validation, and computation of connected components for `TriMesh`.
  - Add the `alloc_instead_of_core`, `std_instead_of_alloc`, and `std_instead_of_core` Clippy lints to the workspace.
  - Use `core` and `alloc` directly rather than using an `std` alias.
  - Use `hashbrown` instead of `rustc-hash` when `enhanced-determinism` is not enabled.
  - Make `spade` optional.

### Fix

- Fix trimesh inertia tensor computation [#331](https://github.com/dimforge/parry/pull/331).
- Fix shifted inertia tensor computation [#334](https://github.com/dimforge/parry/pull/334).

## v0.18.0

### Added

- Implement `::to_trimesh` in 2d for `Cuboid` and `Aabb`.
- Fix some edge-cases in `point_in_poly2d` for self-intersecting polygons.
- Fix some edge-cases in mesh/mesh intersection that could result in degenerate triangles being generated.

### Fix

- Fix panic in `epa3::EPA::closest_points` and `epa2::EPA::closest_points`. Related issues: [#253](https://github.com/dimforge/parry/issues/253), [#246](https://github.com/dimforge/parry/issues/246)

### Modified

- Propagate error information while creating a mesh and using functions making use of it (See #262):
  - `TriMesh::new`
  - `TriMesh::intersection_with_aabb`
  - `SharedShape::trimesh`
  - `SharedShape::trimesh_with_flags`
- `point_cloud_bounding_sphere` and `point_cloud_bounding_sphere_with_center` now returns a `BoundingSphere`.
- Removed `IntersectionCompositeShapeShapeBestFirstVisitor` (which had been deprecated for a while):
  use `IntersectionCompositeShapeShapeVisitor` instead.

## v0.17.5

### Fix

- Always compute connected-components from union-find instead of topology. It is faster and the function based on
  topology could result in a crash for non-manifold meshes.

## v0.17.4

### Added

- Add `TriMeshConnectedComponents::to_meshes` and `::to_mesh_buffers` to easily extract individual meshes from the set
  of connected components.
- Add `TriMesh::connected_component_meshes` to get the connected components as meshes directly.

### Modified

- Connected-components extraction will never fail now, and no longer require the successful calculation of the mesh’s
  half-edge topology.

## v0.17.3

### Fix

- Fix compiling with `enhanced-determinism` feature enabled.
  - This is now checked on CI.

## v0.17.2

### Added

- Implement `Shape::feature_normal_at_point` for `TriMesh` to retrieve the normal of a face, when passing a
  `FeatureId::Face`.
- Add `convex_polygons_intersection_points_with_tolerances`, `convex_polygons_intersection_with_tolerances`, and
  `intersect_meshes_with_tolerances` that let the user specify tolerances value for the collinearity check.

### Fix

- Fix some robustness issues in mesh/mesh intersection when parts of both meshes overlap perfectly.
- Improve robustness of convex polygons intersections when all the vertices of one polygon are located in either the
  edges or vertices of the other polygon.
- Fix incorrect orientation sometimes given to the polygon output by the convex polygon intersections when one of the
  polygon is completely inside the other.

## v0.17.1

### Modified

- Improve convergence of epa algorithm in degenerate configurations.
- Fix bug in the mesh/mesh intersection algorithm that didn’t properly take mesh transforms into account.

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
  `query::time_of_impact`) now take a `ShapeCastOptions` instead of the `max_toi` and
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
  `Qbvh::traverse_bvtt_node_parallel`, `Qbvh::traverse_depth_first_parallel`,
  `Qbvh::traverse_depth_first_node_parallel`.

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
