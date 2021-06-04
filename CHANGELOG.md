# Change Log

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
- Renamed `SimdQuadTree` to `QBVH` (Quaternary Bounding Volume Hierarchy). The
  incorrect name `SimdQuadTree` is now deprecated.
- `QBVH::clear_and_rebuild` is now slightly more general.


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
- Add `AABB::EDGES_VERTEX_IDS` and `AABB::FACES_VERTEX_IDS` which are index tables describing
  the vertices of a given edge, or of a given face, of the AABB.

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
- Add `AABB::vertices` to get an array containing its vertices.
- Add `AABB::split_at_center` to split an AABB into 4 (in 2D) or 8 (in 3D) parts.
- Add `AABB::to_trimesh` to compute the mesh representation of an AABB.

### Removed
- Remove the `Shape::as_serialize` method.

### Fixed
- Fix a bug causing making some ball/convex shape contact manifold computation
  fail when they are penetrating deeply.