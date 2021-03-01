# Change Log
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