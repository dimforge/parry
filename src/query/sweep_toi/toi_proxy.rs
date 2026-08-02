use crate::bounding_volume::Aabb;
use crate::math::{Pose, Real, Vector};
use crate::shape::{Shape, TypedShape};

/// Maximum number of proxy points stored inline (cuboid corners).
#[cfg(feature = "dim2")]
pub const TOI_PROXY_INLINE_POINTS: usize = 4;
/// Maximum number of proxy points stored inline (cuboid corners).
#[cfg(feature = "dim3")]
pub const TOI_PROXY_INLINE_POINTS: usize = 8;

enum ProxyPoints<'a> {
    Inline([Vector; TOI_PROXY_INLINE_POINTS], u8),
    Borrowed(&'a [Vector]),
}

/// A point cloud with a radius, approximating a convex shape for sweep-based time-of-impact
/// computations.
///
/// Only shapes that decompose exactly into a point cloud plus an inflation radius can be
/// represented (balls, capsules, segments, triangles, cuboids, convex polygons/polyhedra and
/// their round variants). Other shapes (cylinders, cones, half-spaces, composites, custom
/// shapes) return `None` from [`ToiProxy::from_shape`].
pub struct ToiProxy<'a> {
    points: ProxyPoints<'a>,
    /// The inflation radius around the point cloud.
    pub radius: Real,
}

impl<'a> ToiProxy<'a> {
    /// A proxy made of a single point with a radius.
    pub fn point(point: Vector, radius: Real) -> Self {
        let mut buf = [Vector::ZERO; TOI_PROXY_INLINE_POINTS];
        buf[0] = point;
        Self {
            points: ProxyPoints::Inline(buf, 1),
            radius,
        }
    }

    /// A proxy borrowing its point cloud.
    pub fn from_points(points: &'a [Vector], radius: Real) -> Self {
        assert!(!points.is_empty());
        Self {
            points: ProxyPoints::Borrowed(points),
            radius,
        }
    }

    /// A proxy from an inline array of points.
    pub fn from_array<const N: usize>(points: [Vector; N], radius: Real) -> Self {
        assert!(N > 0 && N <= TOI_PROXY_INLINE_POINTS);
        let mut buf = [Vector::ZERO; TOI_PROXY_INLINE_POINTS];
        buf[..N].copy_from_slice(&points);
        Self {
            points: ProxyPoints::Inline(buf, N as u8),
            radius,
        }
    }

    /// Extracts a proxy from a shape, if the shape decomposes into points + radius.
    pub fn from_shape(shape: &'a dyn Shape) -> Option<Self> {
        match shape.as_typed_shape() {
            TypedShape::Ball(ball) => Some(Self::point(Vector::ZERO, ball.radius)),
            TypedShape::Cuboid(cuboid) => {
                Some(Self::from_cuboid_half_extents(cuboid.half_extents, 0.0))
            }
            TypedShape::RoundCuboid(round) => Some(Self::from_cuboid_half_extents(
                round.inner_shape.half_extents,
                round.border_radius,
            )),
            TypedShape::Capsule(capsule) => Some(Self::from_array(
                [capsule.segment.a, capsule.segment.b],
                capsule.radius,
            )),
            TypedShape::Segment(segment) => Some(Self::from_array([segment.a, segment.b], 0.0)),
            TypedShape::Triangle(tri) => Some(Self::from_array([tri.a, tri.b, tri.c], 0.0)),
            TypedShape::RoundTriangle(round) => {
                let tri = &round.inner_shape;
                Some(Self::from_array([tri.a, tri.b, tri.c], round.border_radius))
            }
            #[cfg(feature = "dim2")]
            #[cfg(feature = "alloc")]
            TypedShape::ConvexPolygon(poly) => Some(Self::from_points(poly.points(), 0.0)),
            #[cfg(feature = "dim2")]
            #[cfg(feature = "alloc")]
            TypedShape::RoundConvexPolygon(round) => Some(Self::from_points(
                round.inner_shape.points(),
                round.border_radius,
            )),
            #[cfg(feature = "dim3")]
            #[cfg(feature = "alloc")]
            TypedShape::ConvexPolyhedron(poly) => Some(Self::from_points(poly.points(), 0.0)),
            #[cfg(feature = "dim3")]
            #[cfg(feature = "alloc")]
            TypedShape::RoundConvexPolyhedron(round) => Some(Self::from_points(
                round.inner_shape.points(),
                round.border_radius,
            )),
            _ => None,
        }
    }

    fn from_cuboid_half_extents(he: Vector, radius: Real) -> Self {
        #[cfg(feature = "dim2")]
        {
            Self::from_array(
                [
                    Vector::new(-he.x, -he.y),
                    Vector::new(he.x, -he.y),
                    Vector::new(he.x, he.y),
                    Vector::new(-he.x, he.y),
                ],
                radius,
            )
        }
        #[cfg(feature = "dim3")]
        {
            Self::from_array(
                [
                    Vector::new(-he.x, -he.y, -he.z),
                    Vector::new(he.x, -he.y, -he.z),
                    Vector::new(he.x, he.y, -he.z),
                    Vector::new(-he.x, he.y, -he.z),
                    Vector::new(-he.x, -he.y, he.z),
                    Vector::new(he.x, -he.y, he.z),
                    Vector::new(he.x, he.y, he.z),
                    Vector::new(-he.x, he.y, he.z),
                ],
                radius,
            )
        }
    }

    /// The proxy’s point cloud.
    #[inline]
    pub fn points(&self) -> &[Vector] {
        match &self.points {
            ProxyPoints::Inline(buf, len) => &buf[..*len as usize],
            ProxyPoints::Borrowed(points) => points,
        }
    }

    /// The index of the proxy point with the greatest projection along `direction`.
    ///
    /// Projections are measured relative to the first point for better accuracy far from the
    /// origin.
    #[inline]
    pub fn support(&self, direction: Vector) -> u32 {
        let points = self.points();
        let origin = points[0];
        let mut best_index = 0;
        let mut best_value = 0.0;
        for (i, pt) in points.iter().enumerate().skip(1) {
            let value = direction.dot(*pt - origin);
            if value > best_value {
                best_index = i;
                best_value = value;
            }
        }
        best_index as u32
    }

    /// The axis-aligned bounding box of this proxy under the given pose.
    pub fn compute_aabb(&self, pose: &Pose) -> Aabb {
        let points = self.points();
        let mut mins = pose.transform_point(points[0]);
        let mut maxs = mins;
        for pt in &points[1..] {
            let p = pose.transform_point(*pt);
            mins = mins.min(p);
            maxs = maxs.max(p);
        }
        Aabb::new(
            mins - Vector::splat(self.radius),
            maxs + Vector::splat(self.radius),
        )
    }
}
