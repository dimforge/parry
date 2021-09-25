use crate::math::{Isometry, Point, Real, Vector, DIM};
#[cfg(feature = "dim2")]
use crate::shape::ConvexPolygon;
#[cfg(feature = "serde-serialize")]
use crate::shape::DeserializableTypedShape;
use crate::shape::{
    Ball, Capsule, Compound, Cuboid, HalfSpace, HeightField, Polyline, RoundShape, Segment, Shape,
    TriMesh, Triangle,
};
#[cfg(feature = "dim3")]
use crate::shape::{Cone, ConvexPolyhedron, Cylinder};
use crate::transformation::vhacd::{VHACDParameters, VHACD};
use na::Unit;
use std::ops::Deref;
use std::sync::Arc;

/// The shape of a collider.
#[derive(Clone)]
pub struct SharedShape(pub Arc<dyn Shape>);

impl Deref for SharedShape {
    type Target = dyn Shape;
    fn deref(&self) -> &dyn Shape {
        &*self.0
    }
}

impl AsRef<dyn Shape> for SharedShape {
    fn as_ref(&self) -> &dyn Shape {
        &*self.0
    }
}

impl SharedShape {
    /// Wraps the given shape as a shared shape.
    pub fn new(shape: impl Shape) -> Self {
        Self(Arc::new(shape))
    }

    /// If this shape is shared, then the content of `self` is cloned into a unique instance,
    /// and a mutable reference to that instance is returned.
    pub fn make_mut(&mut self) -> &mut dyn Shape {
        if Arc::get_mut(&mut self.0).is_none() {
            let unique_self = self.0.clone_box();
            self.0 = unique_self.into();
        }
        Arc::get_mut(&mut self.0).unwrap()
    }

    /// Initialize a compound shape defined by its subshapes.
    pub fn compound(shapes: Vec<(Isometry<Real>, SharedShape)>) -> Self {
        let raw_shapes = shapes.into_iter().map(|s| (s.0, s.1)).collect();
        let compound = Compound::new(raw_shapes);
        SharedShape(Arc::new(compound))
    }

    /// Initialize a ball shape defined by its radius.
    pub fn ball(radius: Real) -> Self {
        SharedShape(Arc::new(Ball::new(radius)))
    }

    /// Initialize a plane shape defined by its outward normal.
    pub fn halfspace(outward_normal: Unit<Vector<Real>>) -> Self {
        SharedShape(Arc::new(HalfSpace::new(outward_normal)))
    }

    /// Initialize a cylindrical shape defined by its half-height
    /// (along along the y axis) and its radius.
    #[cfg(feature = "dim3")]
    pub fn cylinder(half_height: Real, radius: Real) -> Self {
        SharedShape(Arc::new(Cylinder::new(half_height, radius)))
    }

    /// Initialize a rounded cylindrical shape defined by its half-height
    /// (along along the y axis), its radius, and its roundedness (the
    /// radius of the sphere used for dilating the cylinder).
    #[cfg(feature = "dim3")]
    pub fn round_cylinder(half_height: Real, radius: Real, border_radius: Real) -> Self {
        SharedShape(Arc::new(RoundShape {
            base_shape: Cylinder::new(half_height, radius),
            border_radius,
        }))
    }

    /// Initialize a rounded cone shape defined by its half-height
    /// (along along the y axis), its radius, and its roundedness (the
    /// radius of the sphere used for dilating the cylinder).
    #[cfg(feature = "dim3")]
    pub fn round_cone(half_height: Real, radius: Real, border_radius: Real) -> Self {
        SharedShape(Arc::new(RoundShape {
            base_shape: Cone::new(half_height, radius),
            border_radius,
        }))
    }

    /// Initialize a cone shape defined by its half-height
    /// (along along the y axis) and its basis radius.
    #[cfg(feature = "dim3")]
    pub fn cone(half_height: Real, radius: Real) -> Self {
        SharedShape(Arc::new(Cone::new(half_height, radius)))
    }

    /// Initialize a cuboid shape defined by its half-extents.
    #[cfg(feature = "dim2")]
    pub fn cuboid(hx: Real, hy: Real) -> Self {
        SharedShape(Arc::new(Cuboid::new(Vector::new(hx, hy))))
    }

    /// Initialize a round cuboid shape defined by its half-extents and border radius.
    #[cfg(feature = "dim2")]
    pub fn round_cuboid(hx: Real, hy: Real, border_radius: Real) -> Self {
        SharedShape(Arc::new(RoundShape {
            base_shape: Cuboid::new(Vector::new(hx, hy)),
            border_radius,
        }))
    }

    /// Initialize a cuboid shape defined by its half-extents.
    #[cfg(feature = "dim3")]
    pub fn cuboid(hx: Real, hy: Real, hz: Real) -> Self {
        SharedShape(Arc::new(Cuboid::new(Vector::new(hx, hy, hz))))
    }

    /// Initialize a round cuboid shape defined by its half-extents and border radius.
    #[cfg(feature = "dim3")]
    pub fn round_cuboid(hx: Real, hy: Real, hz: Real, border_radius: Real) -> Self {
        SharedShape(Arc::new(RoundShape {
            base_shape: Cuboid::new(Vector::new(hx, hy, hz)),
            border_radius,
        }))
    }

    /// Initialize a capsule shape from its endpoints and radius.
    pub fn capsule(a: Point<Real>, b: Point<Real>, radius: Real) -> Self {
        SharedShape(Arc::new(Capsule::new(a, b, radius)))
    }

    /// Initialize a segment shape from its endpoints.
    pub fn segment(a: Point<Real>, b: Point<Real>) -> Self {
        SharedShape(Arc::new(Segment::new(a, b)))
    }

    /// Initializes a triangle shape.
    pub fn triangle(a: Point<Real>, b: Point<Real>, c: Point<Real>) -> Self {
        SharedShape(Arc::new(Triangle::new(a, b, c)))
    }
    /// Initializes a triangle shape with round corners.
    pub fn round_triangle(
        a: Point<Real>,
        b: Point<Real>,
        c: Point<Real>,
        border_radius: Real,
    ) -> Self {
        SharedShape(Arc::new(RoundShape {
            base_shape: Triangle::new(a, b, c),
            border_radius,
        }))
    }

    /// Initializes a polyline shape defined by its vertex and index buffers.
    ///
    /// If no index buffer is provided, the polyline is assumed to describe a line strip.
    pub fn polyline(vertices: Vec<Point<Real>>, indices: Option<Vec<[u32; 2]>>) -> Self {
        SharedShape(Arc::new(Polyline::new(vertices, indices)))
    }

    /// Initializes a triangle mesh shape defined by its vertex and index buffers.
    pub fn trimesh(vertices: Vec<Point<Real>>, indices: Vec<[u32; 3]>) -> Self {
        SharedShape(Arc::new(TriMesh::new(vertices, indices)))
    }

    /// Initializes a compound shape obtained from the decomposition of the given trimesh (in 3D) or
    /// polyline (in 2D) into convex parts.
    pub fn convex_decomposition(vertices: &[Point<Real>], indices: &[[u32; DIM]]) -> Self {
        Self::convex_decomposition_with_params(vertices, indices, &VHACDParameters::default())
    }

    /// Initializes a compound shape obtained from the decomposition of the given trimesh (in 3D) or
    /// polyline (in 2D) into convex parts dilated with round corners.
    pub fn round_convex_decomposition(
        vertices: &[Point<Real>],
        indices: &[[u32; DIM]],
        border_radius: Real,
    ) -> Self {
        Self::round_convex_decomposition_with_params(
            vertices,
            indices,
            &VHACDParameters::default(),
            border_radius,
        )
    }

    /// Initializes a compound shape obtained from the decomposition of the given trimesh (in 3D) or
    /// polyline (in 2D) into convex parts.
    pub fn convex_decomposition_with_params(
        vertices: &[Point<Real>],
        indices: &[[u32; DIM]],
        params: &VHACDParameters,
    ) -> Self {
        let mut parts = vec![];
        let decomp = VHACD::decompose(params, &vertices, &indices, true);

        #[cfg(feature = "dim2")]
        for vertices in decomp.compute_exact_convex_hulls(&vertices, &indices) {
            if let Some(convex) = Self::convex_polyline(vertices) {
                parts.push((Isometry::identity(), convex));
            }
        }

        #[cfg(feature = "dim3")]
        for (vertices, indices) in decomp.compute_exact_convex_hulls(&vertices, &indices) {
            if let Some(convex) = Self::convex_mesh(vertices, &indices) {
                parts.push((Isometry::identity(), convex));
            }
        }

        Self::compound(parts)
    }

    /// Initializes a compound shape obtained from the decomposition of the given trimesh (in 3D) or
    /// polyline (in 2D) into convex parts dilated with round corners.
    pub fn round_convex_decomposition_with_params(
        vertices: &[Point<Real>],
        indices: &[[u32; DIM]],
        params: &VHACDParameters,
        border_radius: Real,
    ) -> Self {
        let mut parts = vec![];
        let decomp = VHACD::decompose(params, &vertices, &indices, true);

        #[cfg(feature = "dim2")]
        for vertices in decomp.compute_exact_convex_hulls(&vertices, &indices) {
            if let Some(convex) = Self::round_convex_polyline(vertices, border_radius) {
                parts.push((Isometry::identity(), convex));
            }
        }

        #[cfg(feature = "dim3")]
        for (vertices, indices) in decomp.compute_exact_convex_hulls(&vertices, &indices) {
            if let Some(convex) = Self::round_convex_mesh(vertices, &indices, border_radius) {
                parts.push((Isometry::identity(), convex));
            }
        }

        Self::compound(parts)
    }

    /// Creates a new shared shape that is the convex-hull of the given points.
    pub fn convex_hull(points: &[Point<Real>]) -> Option<Self> {
        #[cfg(feature = "dim2")]
        return ConvexPolygon::from_convex_hull(points).map(|ch| SharedShape(Arc::new(ch)));
        #[cfg(feature = "dim3")]
        return ConvexPolyhedron::from_convex_hull(points).map(|ch| SharedShape(Arc::new(ch)));
    }

    /// Creates a new shared shape that is a convex polygon formed by the
    /// given set of points assumed to form a convex polyline (no convex-hull will be automatically
    /// computed).
    #[cfg(feature = "dim2")]
    pub fn convex_polyline(points: Vec<Point<Real>>) -> Option<Self> {
        ConvexPolygon::from_convex_polyline(points).map(|ch| SharedShape(Arc::new(ch)))
    }

    /// Creates a new shared shape that is a convex polyhedron formed by the
    /// given set of points assumed to form a convex mesh (no convex-hull will be automatically
    /// computed).
    #[cfg(feature = "dim3")]
    pub fn convex_mesh(points: Vec<Point<Real>>, indices: &[[u32; 3]]) -> Option<Self> {
        ConvexPolyhedron::from_convex_mesh(points, indices).map(|ch| SharedShape(Arc::new(ch)))
    }

    /// Creates a new shared shape with rounded corners that is the
    /// convex-hull of the given points, dilated by `border_radius`.
    pub fn round_convex_hull(points: &[Point<Real>], border_radius: Real) -> Option<Self> {
        #[cfg(feature = "dim2")]
        return ConvexPolygon::from_convex_hull(points).map(|ch| {
            SharedShape(Arc::new(RoundShape {
                base_shape: ch,
                border_radius,
            }))
        });
        #[cfg(feature = "dim3")]
        return ConvexPolyhedron::from_convex_hull(points).map(|ch| {
            SharedShape(Arc::new(RoundShape {
                base_shape: ch,
                border_radius,
            }))
        });
    }

    /// Creates a new shared shape with round corners that is a convex polygon formed by the
    /// given set of points assumed to form a convex polyline (no convex-hull will be automatically
    /// computed).
    #[cfg(feature = "dim2")]
    pub fn round_convex_polyline(points: Vec<Point<Real>>, border_radius: Real) -> Option<Self> {
        ConvexPolygon::from_convex_polyline(points).map(|ch| {
            SharedShape(Arc::new(RoundShape {
                base_shape: ch,
                border_radius,
            }))
        })
    }

    /// Creates a new shared shape with round corners that is a convex polyhedron formed by the
    /// given set of points assumed to form a convex mesh (no convex-hull will be automatically
    /// computed).
    #[cfg(feature = "dim3")]
    pub fn round_convex_mesh(
        points: Vec<Point<Real>>,
        indices: &[[u32; 3]],
        border_radius: Real,
    ) -> Option<Self> {
        ConvexPolyhedron::from_convex_mesh(points, indices).map(|ch| {
            SharedShape(Arc::new(RoundShape {
                base_shape: ch,
                border_radius,
            }))
        })
    }

    /// Initializes an heightfield shape defined by its set of height and a scale
    /// factor along each coordinate axis.
    #[cfg(feature = "dim2")]
    pub fn heightfield(heights: na::DVector<Real>, scale: Vector<Real>) -> Self {
        SharedShape(Arc::new(HeightField::new(heights, scale)))
    }

    /// Initializes an heightfield shape on the x-z plane defined by its set of height and a scale
    /// factor along each coordinate axis.
    #[cfg(feature = "dim3")]
    pub fn heightfield(heights: na::DMatrix<Real>, scale: Vector<Real>) -> Self {
        SharedShape(Arc::new(HeightField::new(heights, scale)))
    }
}

#[cfg(feature = "serde-serialize")]
impl serde::Serialize for SharedShape {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        self.0.as_typed_shape().serialize(serializer)
    }
}

#[cfg(feature = "serde-serialize")]
impl<'de> serde::Deserialize<'de> for SharedShape {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        use crate::serde::de::Error;
        DeserializableTypedShape::deserialize(deserializer)?
            .into_shared_shape()
            .ok_or(D::Error::custom("Cannot deserialize custom shape."))
    }
}
