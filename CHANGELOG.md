# Change Log

## v0.2.1 - WIP
### Added
- Add a special case for the triangle/cuboid edge-edge case in the SAT implementation.
- Add contact manifold computation between a convex shape and a HalfSpace.
- Add a special case (instead of GJK) for computing closest points between two cuboids.
- Add a special case (instead of GJK) for computing closest points between a cuboid and a triangle.
- Add a special case (instead of GJK) for computing the distance between a ball and a convex shape.
- Fix the separation computation between a cuboid and a support-map shape in the SAT implementation.

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