//! Shapes supported by parry.

pub use self::ball::Ball;
pub use self::capsule::Capsule;
#[doc(inline)]
pub use self::composite_shape::TypedSimdCompositeShape;
pub use self::cuboid::Cuboid;
pub use self::feature_id::{FeatureId, PackedFeatureId};
pub use self::half_space::HalfSpace;
pub use self::polygonal_feature_map::PolygonalFeatureMap;
pub use self::round_shape::RoundShape;
pub use self::segment::{Segment, SegmentPointLocation};
#[cfg(feature = "serde-serialize")]
pub(crate) use self::shape::DeserializableTypedShape;
#[doc(inline)]
pub use self::shape::{Shape, ShapeType, TypedShape};
#[doc(inline)]
pub use self::support_map::SupportMap;
pub use self::triangle::{Triangle, TriangleOrientation, TrianglePointLocation};

#[cfg(feature = "std")]
pub use self::{
    composite_shape::SimdCompositeShape, compound::Compound, polyline::Polyline,
    shared_shape::SharedShape,
};

#[cfg(feature = "dim2")]
#[cfg(feature = "std")]
pub use self::convex_polygon::ConvexPolygon;
#[cfg(feature = "dim2")]
pub use self::heightfield2::*;
#[cfg(feature = "dim2")]
pub use self::polygonal_feature2d::PolygonalFeature;

#[cfg(feature = "dim3")]
pub use self::cone::Cone;
#[cfg(feature = "dim3")]
#[cfg(feature = "std")]
pub use self::convex_polyhedron::ConvexPolyhedron;
#[cfg(feature = "dim3")]
pub use self::cylinder::Cylinder;
#[cfg(feature = "dim3")]
pub use self::heightfield3::*;
#[cfg(feature = "dim3")]
pub use self::polygonal_feature3d::PolygonalFeature;
#[cfg(feature = "dim3")]
pub use self::tetrahedron::{Tetrahedron, TetrahedronPointLocation};
pub use self::trimesh::*;
pub use self::trimesh_storage::TriMeshStorage;

/// A cylinder dilated by a sphere (so it has round corners).
#[cfg(feature = "dim3")]
pub type RoundCylinder = RoundShape<Cylinder>;
/// A cone dilated by a sphere (so it has round corners).
#[cfg(feature = "dim3")]
pub type RoundCone = RoundShape<Cone>;
/// A cuboid dilated by a sphere (so it has round corners).
pub type RoundCuboid = RoundShape<Cuboid>;
/// A triangle dilated by a sphere (so it has round corners).
pub type RoundTriangle = RoundShape<Triangle>;
/// A convex polyhedron dilated by a sphere (so it has round corners).
#[cfg(feature = "dim3")]
#[cfg(feature = "std")]
pub type RoundConvexPolyhedron = RoundShape<ConvexPolyhedron>;
/// A convex polygon dilated by a sphere (so it has round corners).
#[cfg(feature = "dim2")]
#[cfg(feature = "std")]
pub type RoundConvexPolygon = RoundShape<ConvexPolygon>;

mod ball;
mod capsule;
#[doc(hidden)]
pub mod composite_shape;
#[cfg(feature = "std")]
mod compound;
mod cuboid;
mod half_space;
#[cfg(feature = "std")]
mod polyline;
mod round_shape;
mod segment;
#[doc(hidden)]
pub mod shape;
#[doc(hidden)]
pub mod support_map;
mod triangle;

#[cfg(feature = "dim2")]
#[cfg(feature = "std")]
mod convex_polygon;
#[cfg(feature = "dim2")]
mod heightfield2;

#[cfg(feature = "dim3")]
mod cone;
#[cfg(feature = "dim3")]
#[cfg(feature = "std")]
mod convex_polyhedron;
#[cfg(feature = "dim3")]
mod cylinder;
#[cfg(feature = "dim3")]
mod heightfield3;
#[cfg(feature = "dim3")]
mod polygonal_feature3d;
mod polygonal_feature_map;
#[cfg(feature = "dim3")]
mod tetrahedron;
pub(crate) mod trimesh;
// TODO: move this elsewhere?
mod feature_id;
#[cfg(feature = "dim2")]
mod polygonal_feature2d;
#[cfg(feature = "std")]
mod shared_shape;
mod trimesh_storage;
