use crate::bounding_volume::{BoundingVolume, AABB};
use crate::mass_properties::MassProperties;
use crate::math::{Isometry, Point, Real, Vector};
use crate::query::{PointQuery, RayCast};
use crate::shape::composite_shape::SimdCompositeShape;
#[cfg(feature = "serde-serialize")]
use crate::shape::SharedShape;
use crate::shape::{
    Ball, Capsule, Compound, Cuboid, FeatureId, HalfSpace, HeightField, PolygonalFeatureMap,
    Polyline, RoundCuboid, RoundShape, RoundTriangle, Segment, SupportMap, TriMesh, Triangle,
};
#[cfg(feature = "dim3")]
use crate::shape::{
    Cone, ConvexPolyhedron, Cylinder, RoundCone, RoundConvexPolyhedron, RoundCylinder,
};
#[cfg(feature = "dim2")]
use crate::shape::{ConvexPolygon, RoundConvexPolygon};
use downcast_rs::{impl_downcast, DowncastSync};
use na::Unit;
use num::Zero;
use num_derive::FromPrimitive;

#[derive(Copy, Clone, Debug, FromPrimitive)]
/// Enum representing the type of a shape.
pub enum ShapeType {
    /// A ball shape.
    Ball = 0,
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
    /// A set of segments.
    Polyline,
    /// A shape representing a full half-space.
    HalfSpace,
    /// A heightfield shape.
    HeightField,
    /// A Compound shape.
    Compound,
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
    /// A custom user-defined shape.
    Custom,
}

#[derive(Copy, Clone)]
#[cfg_attr(feature = "serde-serialize", derive(Serialize))]
/// Enum representing the shape with its actual type
pub enum TypedShape<'a> {
    /// A ball shape.
    Ball(&'a Ball),
    /// A cuboid shape.
    Cuboid(&'a Cuboid),
    /// A capsule shape.
    Capsule(&'a Capsule),
    /// A segment shape.
    Segment(&'a Segment),
    /// A triangle shape.
    Triangle(&'a Triangle),
    /// A triangle mesh shape.
    TriMesh(&'a TriMesh),
    /// A set of segments.
    Polyline(&'a Polyline),
    /// A shape representing a full half-space.
    HalfSpace(&'a HalfSpace),
    /// A heightfield shape.
    HeightField(&'a HeightField),
    /// A Compound shape.
    Compound(&'a Compound),
    #[cfg(feature = "dim2")]
    ConvexPolygon(&'a ConvexPolygon),
    #[cfg(feature = "dim3")]
    /// A convex polyhedron.
    ConvexPolyhedron(&'a ConvexPolyhedron),
    #[cfg(feature = "dim3")]
    /// A cylindrical shape.
    Cylinder(&'a Cylinder),
    #[cfg(feature = "dim3")]
    /// A cylindrical shape.
    Cone(&'a Cone),
    // /// A custom shape type.
    // Custom(u8),
    /// A cuboid with rounded corners.
    RoundCuboid(&'a RoundCuboid),
    /// A triangle with rounded corners.
    RoundTriangle(&'a RoundTriangle),
    // /// A triangle-mesh with rounded corners.
    // RoundedTriMesh,
    // /// An heightfield with rounded corners.
    // RoundedHeightField,
    /// A cylinder with rounded corners.
    #[cfg(feature = "dim3")]
    RoundCylinder(&'a RoundCylinder),
    /// A cone with rounded corners.
    #[cfg(feature = "dim3")]
    RoundCone(&'a RoundCone),
    /// A convex polyhedron with rounded corners.
    #[cfg(feature = "dim3")]
    RoundConvexPolyhedron(&'a RoundConvexPolyhedron),
    /// A convex polygon with rounded corners.
    #[cfg(feature = "dim2")]
    RoundConvexPolygon(&'a RoundConvexPolygon),
    /// A custom user-defined shape with a type identified by a number.
    Custom(u32),
}

#[cfg(feature = "serde-serialize")]
#[derive(Deserialize)]
// NOTE: tha this enum MUST match the `TypedShape` enum.
/// Enum representing the shape with its actual type
pub(crate) enum DeserializableTypedShape {
    /// A ball shape.
    Ball(Ball),
    /// A cuboid shape.
    Cuboid(Cuboid),
    /// A capsule shape.
    Capsule(Capsule),
    /// A segment shape.
    Segment(Segment),
    /// A triangle shape.
    Triangle(Triangle),
    /// A triangle mesh shape.
    TriMesh(TriMesh),
    /// A set of segments.
    Polyline(Polyline),
    /// A shape representing a full half-space.
    HalfSpace(HalfSpace),
    /// A heightfield shape.
    HeightField(HeightField),
    /// A Compound shape.
    Compound(Compound),
    #[cfg(feature = "dim2")]
    ConvexPolygon(ConvexPolygon),
    #[cfg(feature = "dim3")]
    /// A convex polyhedron.
    ConvexPolyhedron(ConvexPolyhedron),
    #[cfg(feature = "dim3")]
    /// A cylindrical shape.
    Cylinder(Cylinder),
    #[cfg(feature = "dim3")]
    /// A cylindrical shape.
    Cone(Cone),
    // /// A custom shape type.
    // Custom(u8),
    /// A cuboid with rounded corners.
    RoundCuboid(RoundCuboid),
    /// A triangle with rounded corners.
    RoundTriangle(RoundTriangle),
    // /// A triangle-mesh with rounded corners.
    // RoundedTriMesh,
    // /// An heightfield with rounded corners.
    // RoundedHeightField,
    /// A cylinder with rounded corners.
    #[cfg(feature = "dim3")]
    RoundCylinder(RoundCylinder),
    /// A cone with rounded corners.
    #[cfg(feature = "dim3")]
    RoundCone(RoundCone),
    /// A convex polyhedron with rounded corners.
    #[cfg(feature = "dim3")]
    RoundConvexPolyhedron(RoundConvexPolyhedron),
    /// A convex polygon with rounded corners.
    #[cfg(feature = "dim2")]
    RoundConvexPolygon(RoundConvexPolygon),
    /// A custom user-defined shape identified by a number.
    Custom(u32),
}

#[cfg(feature = "serde-serialize")]
impl DeserializableTypedShape {
    /// Converts `self` to a `SharedShape` if `self` isn't `Custom`.
    pub fn into_shared_shape(self) -> Option<SharedShape> {
        match self {
            DeserializableTypedShape::Ball(s) => Some(SharedShape::new(s)),
            DeserializableTypedShape::Cuboid(s) => Some(SharedShape::new(s)),
            DeserializableTypedShape::Capsule(s) => Some(SharedShape::new(s)),
            DeserializableTypedShape::Segment(s) => Some(SharedShape::new(s)),
            DeserializableTypedShape::Triangle(s) => Some(SharedShape::new(s)),
            DeserializableTypedShape::TriMesh(s) => Some(SharedShape::new(s)),
            DeserializableTypedShape::Polyline(s) => Some(SharedShape::new(s)),
            DeserializableTypedShape::HalfSpace(s) => Some(SharedShape::new(s)),
            DeserializableTypedShape::HeightField(s) => Some(SharedShape::new(s)),
            DeserializableTypedShape::Compound(s) => Some(SharedShape::new(s)),
            #[cfg(feature = "dim2")]
            DeserializableTypedShape::ConvexPolygon(s) => Some(SharedShape::new(s)),
            #[cfg(feature = "dim3")]
            DeserializableTypedShape::ConvexPolyhedron(s) => Some(SharedShape::new(s)),
            #[cfg(feature = "dim3")]
            DeserializableTypedShape::Cylinder(s) => Some(SharedShape::new(s)),
            #[cfg(feature = "dim3")]
            DeserializableTypedShape::Cone(s) => Some(SharedShape::new(s)),
            DeserializableTypedShape::RoundCuboid(s) => Some(SharedShape::new(s)),
            DeserializableTypedShape::RoundTriangle(s) => Some(SharedShape::new(s)),
            #[cfg(feature = "dim3")]
            DeserializableTypedShape::RoundCylinder(s) => Some(SharedShape::new(s)),
            #[cfg(feature = "dim3")]
            DeserializableTypedShape::RoundCone(s) => Some(SharedShape::new(s)),
            #[cfg(feature = "dim3")]
            DeserializableTypedShape::RoundConvexPolyhedron(s) => Some(SharedShape::new(s)),
            #[cfg(feature = "dim2")]
            DeserializableTypedShape::RoundConvexPolygon(s) => Some(SharedShape::new(s)),
            DeserializableTypedShape::Custom(_) => None,
        }
    }
}

/// Trait implemented by shapes usable by Rapier.
pub trait Shape: RayCast + PointQuery + DowncastSync {
    /// Computes the AABB of this shape.
    fn compute_local_aabb(&self) -> AABB;

    /// Computes the AABB of this shape with the given position.
    fn compute_aabb(&self, position: &Isometry<Real>) -> AABB {
        self.compute_local_aabb().transform_by(position)
    }

    /// Compute the mass-properties of this shape given its uniform density.
    fn mass_properties(&self, density: Real) -> MassProperties;

    /// Gets the type tag of this shape.
    fn shape_type(&self) -> ShapeType;

    /// Gets the underlying shape as an enum.
    fn as_typed_shape(&self) -> TypedShape;

    fn ccd_thickness(&self) -> Real;

    /// Is this shape known to be convex?
    ///
    /// If this returns `true` then `self` is known to be convex.
    /// If this returns `false` then it is not known whether or
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
    fn as_polygonal_feature_map(&self) -> Option<(&dyn PolygonalFeatureMap, Real)> {
        None
    }

    // fn as_rounded(&self) -> Option<&Rounded<Box<AnyShape>>> {
    //     None
    // }

    /// The shape's normal at the given point located on a specific feature.
    fn feature_normal_at_point(
        &self,
        _feature: FeatureId,
        _point: &Point<Real>,
    ) -> Option<Unit<Vector<Real>>> {
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

    /// Converts this abstract shape to a compound shape, if it is one.
    pub fn as_compound(&self) -> Option<&Compound> {
        self.downcast_ref()
    }

    /// Converts this abstract shape to a triangle mesh, if it is one.
    pub fn as_trimesh(&self) -> Option<&TriMesh> {
        self.downcast_ref()
    }

    /// Converts this abstract shape to a polyline, if it is one.
    pub fn as_polyline(&self) -> Option<&Polyline> {
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
    fn compute_local_aabb(&self) -> AABB {
        self.local_aabb()
    }

    fn compute_aabb(&self, position: &Isometry<Real>) -> AABB {
        self.aabb(position)
    }

    fn mass_properties(&self, density: Real) -> MassProperties {
        MassProperties::from_ball(density, self.radius)
    }

    fn ccd_thickness(&self) -> Real {
        self.radius
    }

    fn is_convex(&self) -> bool {
        true
    }

    fn shape_type(&self) -> ShapeType {
        ShapeType::Ball
    }

    fn as_typed_shape(&self) -> TypedShape {
        TypedShape::Ball(self)
    }

    fn as_support_map(&self) -> Option<&dyn SupportMap> {
        Some(self as &dyn SupportMap)
    }
}

impl Shape for Cuboid {
    fn compute_local_aabb(&self) -> AABB {
        self.local_aabb()
    }

    fn compute_aabb(&self, position: &Isometry<Real>) -> AABB {
        self.aabb(position)
    }

    fn mass_properties(&self, density: Real) -> MassProperties {
        MassProperties::from_cuboid(density, self.half_extents)
    }

    fn is_convex(&self) -> bool {
        true
    }

    fn shape_type(&self) -> ShapeType {
        ShapeType::Cuboid
    }

    fn as_typed_shape(&self) -> TypedShape {
        TypedShape::Cuboid(self)
    }

    fn ccd_thickness(&self) -> Real {
        self.half_extents.min()
    }

    fn as_support_map(&self) -> Option<&dyn SupportMap> {
        Some(self as &dyn SupportMap)
    }

    fn as_polygonal_feature_map(&self) -> Option<(&dyn PolygonalFeatureMap, Real)> {
        Some((self as &dyn PolygonalFeatureMap, 0.0))
    }
}

impl Shape for Capsule {
    fn compute_local_aabb(&self) -> AABB {
        self.local_aabb()
    }

    fn compute_aabb(&self, position: &Isometry<Real>) -> AABB {
        self.aabb(position)
    }

    fn mass_properties(&self, density: Real) -> MassProperties {
        MassProperties::from_capsule(density, self.segment.a, self.segment.b, self.radius)
    }

    fn is_convex(&self) -> bool {
        true
    }

    fn shape_type(&self) -> ShapeType {
        ShapeType::Capsule
    }

    fn as_typed_shape(&self) -> TypedShape {
        TypedShape::Capsule(self)
    }

    fn ccd_thickness(&self) -> Real {
        self.radius
    }

    fn as_support_map(&self) -> Option<&dyn SupportMap> {
        Some(self as &dyn SupportMap)
    }

    fn as_polygonal_feature_map(&self) -> Option<(&dyn PolygonalFeatureMap, Real)> {
        Some((&self.segment as &dyn PolygonalFeatureMap, self.radius))
    }
}

impl Shape for Triangle {
    fn compute_local_aabb(&self) -> AABB {
        self.local_aabb()
    }

    fn compute_aabb(&self, position: &Isometry<Real>) -> AABB {
        self.aabb(position)
    }

    fn mass_properties(&self, _density: Real) -> MassProperties {
        #[cfg(feature = "dim2")]
        return MassProperties::from_triangle(_density, &self.a, &self.b, &self.c);
        #[cfg(feature = "dim3")]
        return MassProperties::zero();
    }

    fn is_convex(&self) -> bool {
        true
    }

    fn shape_type(&self) -> ShapeType {
        ShapeType::Triangle
    }

    fn as_typed_shape(&self) -> TypedShape {
        TypedShape::Triangle(self)
    }

    fn ccd_thickness(&self) -> Real {
        // TODO: in 2D use the smallest height of the triangle.
        0.0
    }

    fn as_support_map(&self) -> Option<&dyn SupportMap> {
        Some(self as &dyn SupportMap)
    }

    fn as_polygonal_feature_map(&self) -> Option<(&dyn PolygonalFeatureMap, Real)> {
        Some((self as &dyn PolygonalFeatureMap, 0.0))
    }
}

impl Shape for Segment {
    fn compute_local_aabb(&self) -> AABB {
        self.local_aabb()
    }

    fn compute_aabb(&self, position: &Isometry<Real>) -> AABB {
        self.aabb(position)
    }

    fn mass_properties(&self, _density: Real) -> MassProperties {
        MassProperties::zero()
    }

    fn is_convex(&self) -> bool {
        true
    }

    fn ccd_thickness(&self) -> Real {
        0.0
    }

    fn shape_type(&self) -> ShapeType {
        ShapeType::Segment
    }

    fn as_typed_shape(&self) -> TypedShape {
        TypedShape::Segment(self)
    }

    fn as_support_map(&self) -> Option<&dyn SupportMap> {
        Some(self as &dyn SupportMap)
    }

    fn as_polygonal_feature_map(&self) -> Option<(&dyn PolygonalFeatureMap, Real)> {
        Some((self as &dyn PolygonalFeatureMap, 0.0))
    }
}

impl Shape for Compound {
    fn compute_local_aabb(&self) -> AABB {
        *self.local_aabb()
    }

    fn compute_aabb(&self, position: &Isometry<Real>) -> AABB {
        self.local_aabb().transform_by(position)
    }

    fn mass_properties(&self, density: Real) -> MassProperties {
        MassProperties::from_compound(density, self.shapes())
    }

    fn shape_type(&self) -> ShapeType {
        ShapeType::Compound
    }

    fn as_typed_shape(&self) -> TypedShape {
        TypedShape::Compound(self)
    }

    fn ccd_thickness(&self) -> Real {
        self.shapes()
            .iter()
            .fold(Real::MAX, |curr, (_, s)| curr.min(s.ccd_thickness()))
    }

    fn as_composite_shape(&self) -> Option<&dyn SimdCompositeShape> {
        Some(self as &dyn SimdCompositeShape)
    }
}

impl Shape for Polyline {
    fn compute_local_aabb(&self) -> AABB {
        *self.local_aabb()
    }

    fn compute_aabb(&self, position: &Isometry<Real>) -> AABB {
        self.aabb(position)
    }

    fn mass_properties(&self, _density: Real) -> MassProperties {
        MassProperties::zero()
    }

    fn shape_type(&self) -> ShapeType {
        ShapeType::Polyline
    }

    fn as_typed_shape(&self) -> TypedShape {
        TypedShape::Polyline(self)
    }

    fn ccd_thickness(&self) -> Real {
        0.0
    }

    fn as_composite_shape(&self) -> Option<&dyn SimdCompositeShape> {
        Some(self as &dyn SimdCompositeShape)
    }
}

impl Shape for TriMesh {
    fn compute_local_aabb(&self) -> AABB {
        *self.local_aabb()
    }

    fn compute_aabb(&self, position: &Isometry<Real>) -> AABB {
        self.aabb(position)
    }

    fn mass_properties(&self, _density: Real) -> MassProperties {
        #[cfg(feature = "dim2")]
        return MassProperties::from_trimesh(_density, self.vertices(), self.indices());
        #[cfg(feature = "dim3")]
        return MassProperties::zero();
    }

    fn shape_type(&self) -> ShapeType {
        ShapeType::TriMesh
    }

    fn as_typed_shape(&self) -> TypedShape {
        TypedShape::TriMesh(self)
    }

    fn ccd_thickness(&self) -> Real {
        // TODO: in 2D, return the smallest CCD thickness among triangles?
        0.0
    }

    fn as_composite_shape(&self) -> Option<&dyn SimdCompositeShape> {
        Some(self as &dyn SimdCompositeShape)
    }
}

impl Shape for HeightField {
    fn compute_local_aabb(&self) -> AABB {
        self.local_aabb()
    }

    fn compute_aabb(&self, position: &Isometry<Real>) -> AABB {
        self.aabb(position)
    }

    fn mass_properties(&self, _density: Real) -> MassProperties {
        MassProperties::zero()
    }

    fn shape_type(&self) -> ShapeType {
        ShapeType::HeightField
    }

    fn as_typed_shape(&self) -> TypedShape {
        TypedShape::HeightField(self)
    }

    fn ccd_thickness(&self) -> Real {
        0.0
    }
}

#[cfg(feature = "dim2")]
impl Shape for ConvexPolygon {
    fn compute_local_aabb(&self) -> AABB {
        self.local_aabb()
    }

    fn compute_aabb(&self, position: &Isometry<Real>) -> AABB {
        self.aabb(position)
    }

    fn mass_properties(&self, density: Real) -> MassProperties {
        MassProperties::from_convex_polygon(density, &self.points())
    }

    fn is_convex(&self) -> bool {
        true
    }

    fn shape_type(&self) -> ShapeType {
        ShapeType::ConvexPolygon
    }

    fn as_typed_shape(&self) -> TypedShape {
        TypedShape::ConvexPolygon(self)
    }

    fn ccd_thickness(&self) -> Real {
        // TODO: we should use the OBB instead.
        self.compute_local_aabb().half_extents().min()
    }

    fn as_support_map(&self) -> Option<&dyn SupportMap> {
        Some(self as &dyn SupportMap)
    }

    fn as_polygonal_feature_map(&self) -> Option<(&dyn PolygonalFeatureMap, Real)> {
        Some((self as &dyn PolygonalFeatureMap, 0.0))
    }
}

#[cfg(feature = "dim3")]
impl Shape for ConvexPolyhedron {
    fn compute_local_aabb(&self) -> AABB {
        self.local_aabb()
    }

    fn compute_aabb(&self, position: &Isometry<Real>) -> AABB {
        self.aabb(position)
    }

    fn mass_properties(&self, density: Real) -> MassProperties {
        let (vertices, indices) = self.to_trimesh();
        MassProperties::from_convex_polyhedron(density, &vertices, &indices)
    }

    fn is_convex(&self) -> bool {
        true
    }

    fn shape_type(&self) -> ShapeType {
        ShapeType::ConvexPolyhedron
    }

    fn as_typed_shape(&self) -> TypedShape {
        TypedShape::ConvexPolyhedron(self)
    }

    fn ccd_thickness(&self) -> Real {
        // TODO: we should use the OBB instead.
        self.compute_local_aabb().half_extents().min()
    }

    fn as_support_map(&self) -> Option<&dyn SupportMap> {
        Some(self as &dyn SupportMap)
    }

    fn as_polygonal_feature_map(&self) -> Option<(&dyn PolygonalFeatureMap, Real)> {
        Some((self as &dyn PolygonalFeatureMap, 0.0))
    }
}

#[cfg(feature = "dim3")]
impl Shape for Cylinder {
    fn compute_local_aabb(&self) -> AABB {
        self.local_aabb()
    }

    fn compute_aabb(&self, position: &Isometry<Real>) -> AABB {
        self.aabb(position)
    }

    fn mass_properties(&self, density: Real) -> MassProperties {
        MassProperties::from_cylinder(density, self.half_height, self.radius)
    }

    fn is_convex(&self) -> bool {
        true
    }

    fn shape_type(&self) -> ShapeType {
        ShapeType::Cylinder
    }

    fn as_typed_shape(&self) -> TypedShape {
        TypedShape::Cylinder(self)
    }

    fn ccd_thickness(&self) -> Real {
        self.radius
    }

    fn as_support_map(&self) -> Option<&dyn SupportMap> {
        Some(self as &dyn SupportMap)
    }

    fn as_polygonal_feature_map(&self) -> Option<(&dyn PolygonalFeatureMap, Real)> {
        Some((self as &dyn PolygonalFeatureMap, 0.0))
    }
}

#[cfg(feature = "dim3")]
impl Shape for Cone {
    fn compute_local_aabb(&self) -> AABB {
        self.local_aabb()
    }

    fn compute_aabb(&self, position: &Isometry<Real>) -> AABB {
        self.aabb(position)
    }

    fn mass_properties(&self, density: Real) -> MassProperties {
        MassProperties::from_cone(density, self.half_height, self.radius)
    }

    fn is_convex(&self) -> bool {
        true
    }

    fn shape_type(&self) -> ShapeType {
        ShapeType::Cone
    }

    fn as_typed_shape(&self) -> TypedShape {
        TypedShape::Cone(self)
    }

    fn ccd_thickness(&self) -> Real {
        self.radius
    }

    fn as_support_map(&self) -> Option<&dyn SupportMap> {
        Some(self as &dyn SupportMap)
    }

    fn as_polygonal_feature_map(&self) -> Option<(&dyn PolygonalFeatureMap, Real)> {
        Some((self as &dyn PolygonalFeatureMap, 0.0))
    }
}

impl Shape for HalfSpace {
    fn compute_local_aabb(&self) -> AABB {
        self.local_aabb()
    }

    fn compute_aabb(&self, position: &Isometry<Real>) -> AABB {
        self.aabb(position)
    }

    fn is_convex(&self) -> bool {
        true
    }

    fn ccd_thickness(&self) -> Real {
        f32::MAX as Real
    }

    fn mass_properties(&self, _: Real) -> MassProperties {
        MassProperties::zero()
    }

    fn shape_type(&self) -> ShapeType {
        ShapeType::HalfSpace
    }

    fn as_typed_shape(&self) -> TypedShape {
        TypedShape::HalfSpace(self)
    }
}

macro_rules! impl_shape_for_round_shape(
    ($($S: ty, $Tag: ident);*) => {$(
        impl Shape for RoundShape<$S> {
            fn compute_local_aabb(&self) -> AABB {
                self.base_shape.local_aabb().loosened(self.border_radius)
            }

            fn compute_aabb(&self, position: &Isometry<Real>) -> AABB {
                self.base_shape.aabb(position).loosened(self.border_radius)
            }

            fn mass_properties(&self, density: Real) -> MassProperties {
                self.base_shape.mass_properties(density)
            }

            fn is_convex(&self) -> bool {
                self.base_shape.is_convex()
            }

            fn shape_type(&self) -> ShapeType {
                ShapeType::$Tag
            }

            fn as_typed_shape(&self) -> TypedShape {
                TypedShape::$Tag(self)
            }

            fn ccd_thickness(&self) -> Real {
                self.base_shape.ccd_thickness() + self.border_radius
            }

            fn as_support_map(&self) -> Option<&dyn SupportMap> {
                Some(self as &dyn SupportMap)
            }

            fn as_polygonal_feature_map(&self) -> Option<(&dyn PolygonalFeatureMap, Real)> {
                Some((&self.base_shape as &dyn PolygonalFeatureMap, self.border_radius))
            }
        }
    )*}
);

impl_shape_for_round_shape!(
    Cuboid, RoundCuboid;
    Triangle, RoundTriangle
);
#[cfg(feature = "dim2")]
impl_shape_for_round_shape!(ConvexPolygon, RoundConvexPolygon);
#[cfg(feature = "dim3")]
impl_shape_for_round_shape!(
    Cylinder, RoundCylinder;
    Cone, RoundCone;
    ConvexPolyhedron, RoundConvexPolyhedron
);
