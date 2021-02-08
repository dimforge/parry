# Change Log
## v0.2.0 (WIP)
- Remove the `Shape::as_serialize` method.
- Add the `Shape::as_typed_shape` method to convert the shape trait-object
  to an enum.
- Add the `TypedShape` enum with variants containing references to concrete shapes
  (except for user-defined custom shapes).