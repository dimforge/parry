use crate::bounding_volume::{BoundingVolume, AABB};
use crate::math::{Isometry, Point, Vector};
use crate::query::{PointQuery, RayCast};
use crate::shape::composite_shape::SimdCompositeShape;
use crate::shape::{
    Ball, Capsule, Cuboid, FeatureId, HalfSpace, HeightField, MassProperties, PolygonalFeatureMap,
    RoundCuboid, RoundShape, RoundTriangle, Segment, SupportMap, TriMesh, Triangle,
};
#[cfg(feature = "dim3")]
use crate::shape::{
    Cone, ConvexPolyhedron, Cylinder, RoundCone, RoundConvexPolyhedron, RoundCylinder,
};
#[cfg(feature = "dim2")]
use crate::shape::{ConvexPolygon, RoundConvexPolygon};
use downcast_rs::{impl_downcast, DowncastSync};
#[cfg(feature = "serde-serialize")]
use erased_serde::Serialize;
use na::Unit;
use num::Zero;
use num_derive::FromPrimitive;

#[derive(Copy, Clone, Debug, FromPrimitive)]
/// Enum representing the type of a shape.
pub enum ShapeType {
    /// A ball shape.
    Ball = 0,
    /// A convex polygon shape.
    Polygon,
    /// A cuboid shape.
    Cuboid,
    /// A capsule shape.
    Capsule,
    /// A segment shape.
    Segment,
    /// A triangle shape.
    Triangle,
    /// A triangle mesh shape.
    TriMesh,
    /// A shape representing a full half-space.
    HalfSpace,
    /// A heightfield shape.
    HeightField,
    #[cfg(feature = "dim2")]
    ConvexPolygon,
    #[cfg(feature = "dim3")]
    /// A convex polyhedron.
    ConvexPolyhedron,
    #[cfg(feature = "dim3")]
    /// A cylindrical shape.
    Cylinder,
    #[cfg(feature = "dim3")]
    /// A cylindrical shape.
    Cone,
    // /// A custom shape type.
    // Custom(u8),
    /// A cuboid with rounded corners.
    RoundCuboid,
    /// A triangle with rounded corners.
    RoundTriangle,
    // /// A triangle-mesh with rounded corners.
    // RoundedTriMesh,
    // /// An heightfield with rounded corners.
    // RoundedHeightField,
    /// A cylinder with rounded corners.
    #[cfg(feature = "dim3")]
    RoundCylinder,
    /// A cone with rounded corners.
    #[cfg(feature = "dim3")]
    RoundCone,
    /// A convex polyhedron with rounded corners.
    #[cfg(feature = "dim3")]
    RoundConvexPolyhedron,
    /// A convex polygon with rounded corners.
    #[cfg(feature = "dim2")]
    RoundConvexPolygon,
}

/// Trait implemented by shapes usable by Rapier.
pub trait Shape: RayCast + PointQuery + DowncastSync {
    /// Convert this shape as a serializable entity.
    #[cfg(feature = "serde-serialize")]
    fn as_serialize(&self) -> Option<&dyn Serialize> {
        None
    }

    /// Computes the AABB of this shape.
    fn compute_local_aabb(&self) -> AABB;

    /// Computes the AABB of this shape with the given position.
    fn compute_aabb(&self, position: &Isometry<f32>) -> AABB {
        self.compute_local_aabb().transform_by(position)
    }

    /// Compute the mass-properties of this shape given its uniform density.
    fn mass_properties(&self, density: f32) -> MassProperties;

    /// Gets the type tag of this shape.
    fn shape_type(&self) -> ShapeType;

    /// Is this shape known to be convex?
    ///
    /// If this returns `true` then `self` is known to be convex.
    /// If this returns `false` then `self` it is not known whether or
    /// not `self` is convex.
    fn is_convex(&self) -> bool {
        false
    }

    /// Convents this shape into its support mapping, if it has one.
    fn as_support_map(&self) -> Option<&dyn SupportMap> {
        None
    }

    fn as_composite_shape(&self) -> Option<&dyn SimdCompositeShape> {
        None
    }

    /// Converts this shape to a polygonal feature-map, if it is one.
    fn as_polygonal_feature_map(&self) -> Option<(&dyn PolygonalFeatureMap, f32)> {
        None
    }

    // fn as_rounded(&self) -> Option<&Rounded<Box<AnyShape>>> {
    //     None
    // }

    /// The shape's normal at the given point located on a specific feature.
    fn feature_normal_at_point(
        &self,
        _feature: FeatureId,
        _point: &Point<f32>,
    ) -> Option<Unit<Vector<f32>>> {
        None
    }
}

impl_downcast!(sync Shape);

impl dyn Shape {
    /// Converts this abstract shape to the given shape, if it is one.
    pub fn as_shape<T: Shape>(&self) -> Option<&T> {
        self.downcast_ref()
    }

    /// Converts this abstract shape to a ball, if it is one.
    pub fn as_ball(&self) -> Option<&Ball> {
        self.downcast_ref()
    }

    /// Converts this abstract shape to a cuboid, if it is one.
    pub fn as_cuboid(&self) -> Option<&Cuboid> {
        self.downcast_ref()
    }

    /// Converts this abstract shape to a capsule, if it is one.
    pub fn as_capsule(&self) -> Option<&Capsule> {
        self.downcast_ref()
    }

    /// Converts this abstract shape to a triangle, if it is one.
    pub fn as_triangle(&self) -> Option<&Triangle> {
        self.downcast_ref()
    }

    /// Converts this abstract shape to a triangle mesh, if it is one.
    pub fn as_trimesh(&self) -> Option<&TriMesh> {
        self.downcast_ref()
    }

    /// Converts this abstract shape to a heightfield, if it is one.
    pub fn as_heightfield(&self) -> Option<&HeightField> {
        self.downcast_ref()
    }

    /// Converts this abstract shape to a round cuboid, if it is one.
    pub fn as_round_cuboid(&self) -> Option<&RoundCuboid> {
        self.downcast_ref()
    }

    /// Converts this abstract shape to a round triangle, if it is one.
    pub fn as_round_triangle(&self) -> Option<&RoundTriangle> {
        self.downcast_ref()
    }

    #[cfg(feature = "dim2")]
    pub fn as_convex_polygon(&self) -> Option<&ConvexPolygon> {
        self.downcast_ref()
    }

    /// Converts this abstract shape to a round convex polygon, if it is one.
    #[cfg(feature = "dim2")]
    pub fn as_round_convex_polygon(&self) -> Option<&RoundConvexPolygon> {
        self.downcast_ref()
    }

    #[cfg(feature = "dim3")]
    pub fn as_convex_polyhedron(&self) -> Option<&ConvexPolyhedron> {
        self.downcast_ref()
    }

    /// Converts this abstract shape to a cylinder, if it is one.
    #[cfg(feature = "dim3")]
    pub fn as_cylinder(&self) -> Option<&Cylinder> {
        self.downcast_ref()
    }

    /// Converts this abstract shape to a cone, if it is one.
    #[cfg(feature = "dim3")]
    pub fn as_cone(&self) -> Option<&Cone> {
        self.downcast_ref()
    }

    /// Converts this abstract shape to a round cylinder, if it is one.
    #[cfg(feature = "dim3")]
    pub fn as_round_cylinder(&self) -> Option<&RoundCylinder> {
        self.downcast_ref()
    }

    /// Converts this abstract shape to a round cone, if it is one.
    #[cfg(feature = "dim3")]
    pub fn as_round_cone(&self) -> Option<&RoundCone> {
        self.downcast_ref()
    }

    /// Converts this abstract shape to a round convex polyhedron, if it is one.
    #[cfg(feature = "dim3")]
    pub fn as_round_convex_polyhedron(&self) -> Option<&RoundConvexPolyhedron> {
        self.downcast_ref()
    }
}

impl Shape for Ball {
    #[cfg(feature = "serde-serialize")]
    fn as_serialize(&self) -> Option<&dyn Serialize> {
        Some(self as &dyn Serialize)
    }

    fn compute_local_aabb(&self) -> AABB {
        self.local_aabb()
    }

    fn compute_aabb(&self, position: &Isometry<f32>) -> AABB {
        self.aabb(position)
    }

    fn mass_properties(&self, density: f32) -> MassProperties {
        MassProperties::from_ball(density, self.radius)
    }

    fn is_convex(&self) -> bool {
        true
    }

    fn shape_type(&self) -> ShapeType {
        ShapeType::Ball
    }

    fn as_support_map(&self) -> Option<&dyn SupportMap> {
        Some(self as &dyn SupportMap)
    }
}

// impl Shape for Polygon {
//     #[cfg(feature = "serde-serialize")]
//     fn as_serialize(&self) -> Option<&dyn Serialize> {
//         Some(self as &dyn Serialize)
//     }
//
//     fn compute_aabb(&self, position: &Isometry<f32>) -> AABB {
//         self.aabb(position)
//     }
//
//     fn mass_properties(&self, _density: f32) -> MassProperties {
//         unimplemented!()
//     }
//
//     fn shape_type(&self) -> ShapeType {
//         ShapeType::Polygon
//     }
// }

impl Shape for Cuboid {
    #[cfg(feature = "serde-serialize")]
    fn as_serialize(&self) -> Option<&dyn Serialize> {
        Some(self as &dyn Serialize)
    }

    fn compute_local_aabb(&self) -> AABB {
        self.local_aabb()
    }

    fn compute_aabb(&self, position: &Isometry<f32>) -> AABB {
        self.aabb(position)
    }

    fn mass_properties(&self, density: f32) -> MassProperties {
        MassProperties::from_cuboid(density, self.half_extents)
    }

    fn is_convex(&self) -> bool {
        true
    }

    fn shape_type(&self) -> ShapeType {
        ShapeType::Cuboid
    }

    fn as_support_map(&self) -> Option<&dyn SupportMap> {
        Some(self as &dyn SupportMap)
    }

    fn as_polygonal_feature_map(&self) -> Option<(&dyn PolygonalFeatureMap, f32)> {
        Some((self as &dyn PolygonalFeatureMap, 0.0))
    }
}

impl Shape for Capsule {
    #[cfg(feature = "serde-serialize")]
    fn as_serialize(&self) -> Option<&dyn Serialize> {
        Some(self as &dyn Serialize)
    }

    fn compute_local_aabb(&self) -> AABB {
        self.local_aabb()
    }

    fn compute_aabb(&self, position: &Isometry<f32>) -> AABB {
        self.aabb(position)
    }

    fn mass_properties(&self, density: f32) -> MassProperties {
        MassProperties::from_capsule(density, self.segment.a, self.segment.b, self.radius)
    }

    fn is_convex(&self) -> bool {
        true
    }

    fn shape_type(&self) -> ShapeType {
        ShapeType::Capsule
    }

    fn as_support_map(&self) -> Option<&dyn SupportMap> {
        Some(self as &dyn SupportMap)
    }

    fn as_polygonal_feature_map(&self) -> Option<(&dyn PolygonalFeatureMap, f32)> {
        Some((&self.segment as &dyn PolygonalFeatureMap, self.radius))
    }
}

impl Shape for Triangle {
    #[cfg(feature = "serde-serialize")]
    fn as_serialize(&self) -> Option<&dyn Serialize> {
        Some(self as &dyn Serialize)
    }

    fn compute_local_aabb(&self) -> AABB {
        self.local_aabb()
    }

    fn compute_aabb(&self, position: &Isometry<f32>) -> AABB {
        self.aabb(position)
    }

    fn mass_properties(&self, _density: f32) -> MassProperties {
        MassProperties::zero()
    }

    fn is_convex(&self) -> bool {
        true
    }

    fn shape_type(&self) -> ShapeType {
        ShapeType::Triangle
    }

    fn as_support_map(&self) -> Option<&dyn SupportMap> {
        Some(self as &dyn SupportMap)
    }

    fn as_polygonal_feature_map(&self) -> Option<(&dyn PolygonalFeatureMap, f32)> {
        Some((self as &dyn PolygonalFeatureMap, 0.0))
    }
}

impl Shape for Segment {
    #[cfg(feature = "serde-serialize")]
    fn as_serialize(&self) -> Option<&dyn Serialize> {
        Some(self as &dyn Serialize)
    }

    fn compute_local_aabb(&self) -> AABB {
        self.local_aabb()
    }

    fn compute_aabb(&self, position: &Isometry<f32>) -> AABB {
        self.aabb(position)
    }

    fn mass_properties(&self, _density: f32) -> MassProperties {
        MassProperties::zero()
    }

    fn is_convex(&self) -> bool {
        true
    }

    fn shape_type(&self) -> ShapeType {
        ShapeType::Segment
    }

    fn as_support_map(&self) -> Option<&dyn SupportMap> {
        Some(self as &dyn SupportMap)
    }

    fn as_polygonal_feature_map(&self) -> Option<(&dyn PolygonalFeatureMap, f32)> {
        Some((self as &dyn PolygonalFeatureMap, 0.0))
    }
}

impl Shape for TriMesh {
    #[cfg(feature = "serde-serialize")]
    fn as_serialize(&self) -> Option<&dyn Serialize> {
        Some(self as &dyn Serialize)
    }

    fn compute_local_aabb(&self) -> AABB {
        *self.local_aabb()
    }

    fn compute_aabb(&self, position: &Isometry<f32>) -> AABB {
        self.aabb(position)
    }

    fn mass_properties(&self, _density: f32) -> MassProperties {
        MassProperties::zero()
    }

    fn shape_type(&self) -> ShapeType {
        ShapeType::TriMesh
    }

    fn as_composite_shape(&self) -> Option<&dyn SimdCompositeShape> {
        Some(self as &dyn SimdCompositeShape)
    }
}

impl Shape for HeightField {
    #[cfg(feature = "serde-serialize")]
    fn as_serialize(&self) -> Option<&dyn Serialize> {
        Some(self as &dyn Serialize)
    }

    fn compute_local_aabb(&self) -> AABB {
        self.local_aabb()
    }

    fn compute_aabb(&self, position: &Isometry<f32>) -> AABB {
        self.aabb(position)
    }

    fn mass_properties(&self, _density: f32) -> MassProperties {
        MassProperties::zero()
    }

    fn shape_type(&self) -> ShapeType {
        ShapeType::HeightField
    }
}

#[cfg(feature = "dim2")]
impl Shape for ConvexPolygon {
    #[cfg(feature = "serde-serialize")]
    fn as_serialize(&self) -> Option<&dyn Serialize> {
        Some(self as &dyn Serialize)
    }

    fn compute_local_aabb(&self) -> AABB {
        self.local_aabb()
    }

    fn compute_aabb(&self, position: &Isometry<f32>) -> AABB {
        self.aabb(position)
    }

    fn mass_properties(&self, density: f32) -> MassProperties {
        MassProperties::from_convex_polygon(density, &self.points())
    }

    fn is_convex(&self) -> bool {
        true
    }

    fn shape_type(&self) -> ShapeType {
        ShapeType::ConvexPolygon
    }

    fn as_support_map(&self) -> Option<&dyn SupportMap> {
        Some(self as &dyn SupportMap)
    }

    fn as_polygonal_feature_map(&self) -> Option<(&dyn PolygonalFeatureMap, f32)> {
        Some((self as &dyn PolygonalFeatureMap, 0.0))
    }
}

#[cfg(feature = "dim3")]
impl Shape for ConvexPolyhedron {
    #[cfg(feature = "serde-serialize")]
    fn as_serialize(&self) -> Option<&dyn Serialize> {
        Some(self as &dyn Serialize)
    }

    fn compute_local_aabb(&self) -> AABB {
        self.local_aabb()
    }

    fn compute_aabb(&self, position: &Isometry<f32>) -> AABB {
        self.aabb(position)
    }

    fn mass_properties(&self, density: f32) -> MassProperties {
        let (vertices, indices) = self.to_trimesh();
        MassProperties::from_convex_polyhedron(density, &vertices, &indices)
    }

    fn is_convex(&self) -> bool {
        true
    }

    fn shape_type(&self) -> ShapeType {
        ShapeType::ConvexPolyhedron
    }

    fn as_support_map(&self) -> Option<&dyn SupportMap> {
        Some(self as &dyn SupportMap)
    }

    fn as_polygonal_feature_map(&self) -> Option<(&dyn PolygonalFeatureMap, f32)> {
        Some((self as &dyn PolygonalFeatureMap, 0.0))
    }
}

#[cfg(feature = "dim3")]
impl Shape for Cylinder {
    #[cfg(feature = "serde-serialize")]
    fn as_serialize(&self) -> Option<&dyn Serialize> {
        Some(self as &dyn Serialize)
    }

    fn compute_local_aabb(&self) -> AABB {
        self.local_aabb()
    }

    fn compute_aabb(&self, position: &Isometry<f32>) -> AABB {
        self.aabb(position)
    }

    fn mass_properties(&self, density: f32) -> MassProperties {
        MassProperties::from_cylinder(density, self.half_height, self.radius)
    }

    fn is_convex(&self) -> bool {
        true
    }

    fn shape_type(&self) -> ShapeType {
        ShapeType::Cylinder
    }

    fn as_support_map(&self) -> Option<&dyn SupportMap> {
        Some(self as &dyn SupportMap)
    }

    fn as_polygonal_feature_map(&self) -> Option<(&dyn PolygonalFeatureMap, f32)> {
        Some((self as &dyn PolygonalFeatureMap, 0.0))
    }
}

#[cfg(feature = "dim3")]
impl Shape for Cone {
    #[cfg(feature = "serde-serialize")]
    fn as_serialize(&self) -> Option<&dyn Serialize> {
        Some(self as &dyn Serialize)
    }

    fn compute_local_aabb(&self) -> AABB {
        self.local_aabb()
    }

    fn compute_aabb(&self, position: &Isometry<f32>) -> AABB {
        self.aabb(position)
    }

    fn mass_properties(&self, density: f32) -> MassProperties {
        MassProperties::from_cone(density, self.half_height, self.radius)
    }

    fn is_convex(&self) -> bool {
        true
    }

    fn shape_type(&self) -> ShapeType {
        ShapeType::Cone
    }

    fn as_support_map(&self) -> Option<&dyn SupportMap> {
        Some(self as &dyn SupportMap)
    }

    fn as_polygonal_feature_map(&self) -> Option<(&dyn PolygonalFeatureMap, f32)> {
        Some((self as &dyn PolygonalFeatureMap, 0.0))
    }
}

impl Shape for HalfSpace {
    #[cfg(feature = "serde-serialize")]
    fn as_serialize(&self) -> Option<&dyn Serialize> {
        Some(self as &dyn Serialize)
    }

    fn compute_local_aabb(&self) -> AABB {
        self.local_aabb()
    }

    fn compute_aabb(&self, position: &Isometry<f32>) -> AABB {
        self.aabb(position)
    }

    fn is_convex(&self) -> bool {
        true
    }

    fn mass_properties(&self, _: f32) -> MassProperties {
        MassProperties::zero()
    }

    fn shape_type(&self) -> ShapeType {
        ShapeType::HalfSpace
    }
}

macro_rules! impl_shape_for_shape_with_border(
    ($($S: ty, $Tag: expr);*) => {$(
        impl Shape for RoundShape<$S> {
            #[cfg(feature = "serde-serialize")]
            fn as_serialize(&self) -> Option<&dyn Serialize> {
                Some(self as &dyn Serialize)
            }

            fn compute_local_aabb(&self) -> AABB {
                self.base_shape.local_aabb().loosened(self.border_radius)
            }

            fn compute_aabb(&self, position: &Isometry<f32>) -> AABB {
                self.base_shape.aabb(position).loosened(self.border_radius)
            }

            fn mass_properties(&self, density: f32) -> MassProperties {
                self.base_shape.mass_properties(density)
            }

            fn is_convex(&self) -> bool {
                self.base_shape.is_convex()
            }

            fn shape_type(&self) -> ShapeType {
                $Tag
            }

            fn as_support_map(&self) -> Option<&dyn SupportMap> {
                Some(self as &dyn SupportMap)
            }

            fn as_polygonal_feature_map(&self) -> Option<(&dyn PolygonalFeatureMap, f32)> {
                Some((&self.base_shape as &dyn PolygonalFeatureMap, self.border_radius))
            }
        }
    )*}
);

impl_shape_for_shape_with_border!(
    Cuboid, ShapeType::RoundCuboid;
    Triangle, ShapeType::RoundTriangle
);
#[cfg(feature = "dim2")]
impl_shape_for_shape_with_border!(ConvexPolygon, ShapeType::RoundConvexPolygon);
#[cfg(feature = "dim3")]
impl_shape_for_shape_with_border!(
    Cylinder, ShapeType::RoundCylinder;
    Cone, ShapeType::RoundCone;
    ConvexPolyhedron, ShapeType::RoundConvexPolyhedron
);
