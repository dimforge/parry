//! Shapes supported by cdl.

pub use self::ball::Ball;
pub use self::capsule::Capsule;
#[doc(inline)]
pub use self::composite_shape::SimdCompositeShape;
pub use self::convex_polyhedron::{ConvexPolyhedron, FeatureId};
pub use self::cuboid::Cuboid;
pub use self::half_space::HalfSpace;
pub use self::mass_properties::MassProperties;
pub use self::polyline::Polyline;
pub use self::segment::{Segment, SegmentPointLocation};
#[doc(inline)]
pub use self::shape::{Shape, ShapeType};
#[doc(inline)]
pub use self::support_map::SupportMap;
pub use self::triangle::{Triangle, TrianglePointLocation};

#[cfg(feature = "dim2")]
pub use self::convex_polygon::ConvexPolygon;
#[cfg(feature = "dim2")]
pub use self::convex_polygonal_feature2::ConvexPolygonalFeature;
#[cfg(feature = "dim2")]
pub use self::cuboid_feature2d::{CuboidFeature, CuboidFeatureFace, CuboidFeatureVertex};
#[cfg(feature = "dim2")]
pub use self::heightfield2::HeightField;

#[cfg(feature = "dim3")]
pub use self::cone::Cone;
#[cfg(feature = "dim3")]
pub use self::convex::ConvexHull;
#[cfg(feature = "dim3")]
pub use self::cylinder::Cylinder;
#[cfg(feature = "dim3")]
pub use self::round_cylinder::RoundCylinder;
//#[cfg(feature = "dim3")]
//pub use self::deformable_trimesh::DeformableTriMesh;
#[cfg(feature = "dim3")]
pub use self::convex_polygonal_feature3::{ClippingCache, ConvexPolygonalFeature};
#[cfg(feature = "dim3")]
pub use self::cuboid_feature3d::{
    CuboidFeature, CuboidFeatureEdge, CuboidFeatureFace, CuboidFeatureVertex,
};
#[cfg(feature = "dim3")]
pub use self::heightfield3::{HeightField, HeightFieldCellStatus};
#[cfg(feature = "dim3")]
pub use self::polygonal_feature_map::PolygonalFeatureMap;
#[cfg(feature = "dim3")]
pub use self::polyhedron_feature3d::PolygonalFeature;
#[cfg(feature = "dim3")]
pub use self::tetrahedron::{Tetrahedron, TetrahedronPointLocation};
pub use self::trimesh::TriMesh;

mod ball;
mod capsule;
#[doc(hidden)]
pub mod composite_shape;
mod convex_polyhedron;
mod cuboid;
mod half_space;
mod polyline;
mod segment;
#[doc(hidden)]
pub mod shape;
#[doc(hidden)]
pub mod support_map;
mod triangle;

#[cfg(feature = "dim2")]
mod convex_polygon;
#[cfg(feature = "dim2")]
mod convex_polygonal_feature2;
#[cfg(feature = "dim2")]
mod cuboid_feature2d;
#[cfg(feature = "dim2")]
mod heightfield2;

#[cfg(feature = "dim3")]
mod cone;
#[cfg(feature = "dim3")]
mod convex;
#[cfg(feature = "dim3")]
mod convex_polygonal_feature3;
#[cfg(feature = "dim3")]
mod cuboid_feature3d;
#[cfg(feature = "dim3")]
mod cylinder;
#[cfg(feature = "dim3")]
mod heightfield3;
#[cfg(feature = "dim3")]
mod polygonal_feature_map;
#[cfg(feature = "dim3")]
mod polyhedron_feature3d;
#[cfg(feature = "dim3")]
mod round_cylinder;
#[cfg(feature = "dim3")]
mod tetrahedron;
mod trimesh;
// TODO: move this elsewhere?
mod mass_properties;
