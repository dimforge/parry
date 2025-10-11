//! Axis Aligned Bounding Box.

use crate::bounding_volume::{BoundingSphere, BoundingVolume};
use crate::math::{Isometry, Point, Real, UnitVector, Vector, DIM, TWO_DIM};
use crate::shape::{Cuboid, SupportMap};
use crate::utils::IsometryOps;
use arrayvec::ArrayVec;
use na;
use num::Bounded;

#[cfg(all(feature = "dim3", not(feature = "std")))]
use na::ComplexField;
// for .sin_cos()

use crate::query::{Ray, RayCast};
#[cfg(feature = "rkyv")]
use rkyv::{bytecheck, CheckBytes};

/// An Axis-Aligned Bounding Box (AABB).
///
/// An AABB is the simplest bounding volume, defined by its minimum and maximum corners.
/// It's called "axis-aligned" because its edges are always parallel to the coordinate axes
/// (X, Y, and Z in 3D), making it very fast to test and compute.
///
/// # Structure
///
/// - **mins**: The point with the smallest coordinates on each axis (bottom-left-back corner)
/// - **maxs**: The point with the largest coordinates on each axis (top-right-front corner)
/// - **Invariant**: `mins.x ≤ maxs.x`, `mins.y ≤ maxs.y` (and `mins.z ≤ maxs.z` in 3D)
///
/// # Properties
///
/// - **Axis-aligned**: Edges always parallel to coordinate axes
/// - **Conservative**: May be larger than the actual shape for rotated objects
/// - **Fast**: Intersection tests are very cheap (just coordinate comparisons)
/// - **Hierarchical**: Perfect for spatial data structures (BVH, quadtree, octree)
///
/// # Use Cases
///
/// AABBs are fundamental to collision detection and are used for:
///
/// - **Broad-phase collision detection**: Quickly eliminate distant object pairs
/// - **Spatial partitioning**: Building BVHs, quadtrees, and octrees
/// - **View frustum culling**: Determining what's visible
/// - **Ray tracing acceleration**: Quickly rejecting non-intersecting rays
/// - **Bounding volume for any shape**: Every shape can compute its AABB
///
/// # Performance
///
/// AABBs are the fastest bounding volume for:
/// - Intersection tests: O(1) with just 6 comparisons (3D)
/// - Merging: O(1) with component-wise min/max
/// - Contains test: O(1) with coordinate comparisons
///
/// # Limitations
///
/// - **Rotation invariance**: Must be recomputed when objects rotate
/// - **Tightness**: May waste space for rotated or complex shapes
/// - **No orientation**: Cannot represent oriented bounding boxes (OBB)
///
/// # Example
///
/// ```rust
/// # #[cfg(all(feature = "dim3", feature = "f32"))] {
/// use parry3d::bounding_volume::Aabb;
/// use nalgebra::Point3;
///
/// // Create an AABB for a unit cube centered at origin
/// let mins = Point3::new(-0.5, -0.5, -0.5);
/// let maxs = Point3::new(0.5, 0.5, 0.5);
/// let aabb = Aabb::new(mins, maxs);
///
/// // Check if a point is inside
/// let point = Point3::new(0.0, 0.0, 0.0);
/// assert!(aabb.contains_local_point(&point));
///
/// // Get center and extents
/// assert_eq!(aabb.center(), Point3::origin());
/// assert_eq!(aabb.extents().x, 1.0); // Full width
/// assert_eq!(aabb.half_extents().x, 0.5); // Half width
/// # }
/// ```
///
/// ```rust
/// # #[cfg(all(feature = "dim3", feature = "f32"))] {
/// use parry3d::bounding_volume::Aabb;
/// use nalgebra::Point3;
///
/// // Create from a set of points
/// let points = vec![
///     Point3::new(1.0, 2.0, 3.0),
///     Point3::new(-1.0, 4.0, 2.0),
///     Point3::new(0.0, 0.0, 5.0),
/// ];
/// let aabb = Aabb::from_points(points);
///
/// // The AABB encloses all points
/// assert_eq!(aabb.mins, Point3::new(-1.0, 0.0, 2.0));
/// assert_eq!(aabb.maxs, Point3::new(1.0, 4.0, 5.0));
/// # }
/// ```
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[cfg_attr(feature = "bytemuck", derive(bytemuck::Pod, bytemuck::Zeroable))]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Deserialize, rkyv::Serialize, CheckBytes),
    archive(as = "Self")
)]
#[derive(Debug, PartialEq, Copy, Clone)]
#[repr(C)]
pub struct Aabb {
    /// The point with minimum coordinates (bottom-left-back corner).
    ///
    /// Each component (`x`, `y`, `z`) should be less than or equal to the
    /// corresponding component in `maxs`.
    pub mins: Point<Real>,

    /// The point with maximum coordinates (top-right-front corner).
    ///
    /// Each component (`x`, `y`, `z`) should be greater than or equal to the
    /// corresponding component in `mins`.
    pub maxs: Point<Real>,
}

impl Aabb {
    /// The vertex indices of each edge of this `Aabb`.
    ///
    /// This gives, for each edge of this `Aabb`, the indices of its
    /// vertices when taken from the `self.vertices()` array.
    /// Here is how the faces are numbered, assuming
    /// a right-handed coordinate system:
    ///
    /// ```text
    ///    y             3 - 2
    ///    |           7 − 6 |
    ///    ___ x       |   | 1  (the zero is below 3 and on the left of 1,
    ///   /            4 - 5     hidden by the 4-5-6-7 face.)
    ///  z
    /// ```
    #[cfg(feature = "dim3")]
    pub const EDGES_VERTEX_IDS: [(usize, usize); 12] = [
        (0, 1),
        (1, 2),
        (3, 2),
        (0, 3),
        (4, 5),
        (5, 6),
        (7, 6),
        (4, 7),
        (0, 4),
        (1, 5),
        (2, 6),
        (3, 7),
    ];

    /// The vertex indices of each face of this `Aabb`.
    ///
    /// This gives, for each face of this `Aabb`, the indices of its
    /// vertices when taken from the `self.vertices()` array.
    /// Here is how the faces are numbered, assuming
    /// a right-handed coordinate system:
    ///
    /// ```text
    ///    y             3 - 2
    ///    |           7 − 6 |
    ///    ___ x       |   | 1  (the zero is below 3 and on the left of 1,
    ///   /            4 - 5     hidden by the 4-5-6-7 face.)
    ///  z
    /// ```
    #[cfg(feature = "dim3")]
    pub const FACES_VERTEX_IDS: [(usize, usize, usize, usize); 6] = [
        // Face with normal +X
        (1, 2, 6, 5),
        // Face with normal -X
        (0, 3, 7, 4),
        // Face with normal +Y
        (2, 3, 7, 6),
        // Face with normal -Y
        (1, 0, 4, 5),
        // Face with normal +Z
        (4, 5, 6, 7),
        // Face with normal -Z
        (0, 1, 2, 3),
    ];

    /// The vertex indices of each face of this `Aabb`.
    ///
    /// This gives, for each face of this `Aabb`, the indices of its
    /// vertices when taken from the `self.vertices()` array.
    /// Here is how the faces are numbered, assuming
    /// a right-handed coordinate system:
    ///
    /// ```text
    ///    y             3 - 2
    ///    |             |   |
    ///    ___ x         0 - 1
    /// ```
    #[cfg(feature = "dim2")]
    pub const FACES_VERTEX_IDS: [(usize, usize); 4] = [
        // Face with normal +X
        (1, 2),
        // Face with normal -X
        (3, 0),
        // Face with normal +Y
        (2, 3),
        // Face with normal -Y
        (0, 1),
    ];

    /// Creates a new AABB from its minimum and maximum corners.
    ///
    /// # Arguments
    ///
    /// * `mins` - The point with the smallest coordinates on each axis
    /// * `maxs` - The point with the largest coordinates on each axis
    ///
    /// # Invariant
    ///
    /// Each component of `mins` should be ≤ the corresponding component of `maxs`.
    ///
    /// # Example
    ///
    /// ```rust
    /// # #[cfg(all(feature = "dim3", feature = "f32"))] {
    /// use parry3d::bounding_volume::Aabb;
    /// use nalgebra::Point3;
    ///
    /// // Create a 2x2x2 cube centered at origin
    /// let aabb = Aabb::new(
    ///     Point3::new(-1.0, -1.0, -1.0),
    ///     Point3::new(1.0, 1.0, 1.0)
    /// );
    ///
    /// assert_eq!(aabb.center(), Point3::origin());
    /// assert_eq!(aabb.extents(), nalgebra::Vector3::new(2.0, 2.0, 2.0));
    /// # }
    /// ```
    #[inline]
    pub fn new(mins: Point<Real>, maxs: Point<Real>) -> Aabb {
        Aabb { mins, maxs }
    }

    /// Creates an invalid AABB with inverted bounds.
    ///
    /// The resulting AABB has `mins` set to maximum values and `maxs` set to
    /// minimum values. This is useful as an initial value for AABB merging
    /// algorithms (similar to starting a min operation with infinity).
    ///
    /// # Example
    ///
    /// ```rust
    /// # #[cfg(all(feature = "dim3", feature = "f32"))]
    /// use parry3d::bounding_volume::{Aabb, BoundingVolume};
    /// use nalgebra::Point3;
    ///
    /// let mut aabb = Aabb::new_invalid();
    ///
    /// // Merge with actual points to build proper AABB
    /// aabb.merge(&Aabb::new(Point3::new(1.0, 2.0, 3.0), Point3::new(1.0, 2.0, 3.0)));
    /// aabb.merge(&Aabb::new(Point3::new(-1.0, 0.0, 2.0), Point3::new(-1.0, 0.0, 2.0)));
    ///
    /// // Now contains the merged bounds
    /// assert_eq!(aabb.mins, Point3::new(-1.0, 0.0, 2.0));
    /// assert_eq!(aabb.maxs, Point3::new(1.0, 2.0, 3.0));
    /// # }
    /// ```
    #[inline]
    pub fn new_invalid() -> Self {
        Self::new(
            Vector::repeat(Real::max_value()).into(),
            Vector::repeat(-Real::max_value()).into(),
        )
    }

    /// Creates a new AABB from its center and half-extents.
    ///
    /// This is often more intuitive than specifying min and max corners.
    ///
    /// # Arguments
    ///
    /// * `center` - The center point of the AABB
    /// * `half_extents` - Half the dimensions along each axis
    ///
    /// # Example
    ///
    /// ```rust
    /// # #[cfg(all(feature = "dim3", feature = "f32"))] {
    /// use parry3d::bounding_volume::Aabb;
    /// use nalgebra::{Point3, Vector3};
    ///
    /// // Create a 10x6x8 box centered at (5, 0, 0)
    /// let aabb = Aabb::from_half_extents(
    ///     Point3::new(5.0, 0.0, 0.0),
    ///     Vector3::new(5.0, 3.0, 4.0)
    /// );
    ///
    /// assert_eq!(aabb.mins, Point3::new(0.0, -3.0, -4.0));
    /// assert_eq!(aabb.maxs, Point3::new(10.0, 3.0, 4.0));
    /// assert_eq!(aabb.center(), Point3::new(5.0, 0.0, 0.0));
    /// # }
    /// ```
    #[inline]
    pub fn from_half_extents(center: Point<Real>, half_extents: Vector<Real>) -> Self {
        Self::new(center - half_extents, center + half_extents)
    }

    /// Creates a new AABB that tightly encloses a set of points (references).
    ///
    /// Computes the minimum and maximum coordinates across all points.
    ///
    /// # Arguments
    ///
    /// * `pts` - An iterator over point references
    ///
    /// # Example
    ///
    /// ```rust
    /// # #[cfg(all(feature = "dim3", feature = "f32"))] {
    /// use parry3d::bounding_volume::Aabb;
    /// use nalgebra::Point3;
    ///
    /// let points = vec![
    ///     Point3::new(1.0, 2.0, 3.0),
    ///     Point3::new(-1.0, 4.0, 2.0),
    ///     Point3::new(0.0, 0.0, 5.0),
    /// ];
    ///
    /// let aabb = Aabb::from_points_ref(&points);
    /// assert_eq!(aabb.mins, Point3::new(-1.0, 0.0, 2.0));
    /// assert_eq!(aabb.maxs, Point3::new(1.0, 4.0, 5.0));
    /// # }
    /// ```
    pub fn from_points_ref<'a, I>(pts: I) -> Self
    where
        I: IntoIterator<Item = &'a Point<Real>>,
    {
        super::aabb_utils::local_point_cloud_aabb(pts.into_iter().copied())
    }

    /// Creates a new AABB that tightly encloses a set of points (values).
    ///
    /// Computes the minimum and maximum coordinates across all points.
    ///
    /// # Arguments
    ///
    /// * `pts` - An iterator over point values
    ///
    /// # Example
    ///
    /// ```rust
    /// # #[cfg(all(feature = "dim3", feature = "f32"))] {
    /// use parry3d::bounding_volume::Aabb;
    /// use nalgebra::Point3;
    ///
    /// let aabb = Aabb::from_points(vec![
    ///     Point3::new(1.0, 2.0, 3.0),
    ///     Point3::new(-1.0, 4.0, 2.0),
    ///     Point3::new(0.0, 0.0, 5.0),
    /// ]);
    ///
    /// assert_eq!(aabb.mins, Point3::new(-1.0, 0.0, 2.0));
    /// assert_eq!(aabb.maxs, Point3::new(1.0, 4.0, 5.0));
    /// # }
    /// ```
    pub fn from_points<I>(pts: I) -> Self
    where
        I: IntoIterator<Item = Point<Real>>,
    {
        super::aabb_utils::local_point_cloud_aabb(pts)
    }

    /// Returns the center point of this AABB.
    ///
    /// The center is the midpoint between `mins` and `maxs`.
    ///
    /// # Example
    ///
    /// ```rust
    /// # #[cfg(all(feature = "dim3", feature = "f32"))] {
    /// use parry3d::bounding_volume::Aabb;
    /// use nalgebra::Point3;
    ///
    /// let aabb = Aabb::new(
    ///     Point3::new(-2.0, -3.0, -4.0),
    ///     Point3::new(2.0, 3.0, 4.0)
    /// );
    ///
    /// assert_eq!(aabb.center(), Point3::origin());
    /// # }
    /// ```
    #[inline]
    pub fn center(&self) -> Point<Real> {
        na::center(&self.mins, &self.maxs)
    }

    /// Returns the half-extents of this AABB.
    ///
    /// Half-extents are half the dimensions along each axis.
    ///
    /// # Example
    ///
    /// ```rust
    /// # #[cfg(all(feature = "dim3", feature = "f32"))] {
    /// use parry3d::bounding_volume::Aabb;
    /// use nalgebra::{Point3, Vector3};
    ///
    /// let aabb = Aabb::new(
    ///     Point3::new(-5.0, -3.0, -2.0),
    ///     Point3::new(5.0, 3.0, 2.0)
    /// );
    ///
    /// let half = aabb.half_extents();
    /// assert_eq!(half, Vector3::new(5.0, 3.0, 2.0));
    ///
    /// // Full dimensions are 2 * half_extents
    /// let full = aabb.extents();
    /// assert_eq!(full, Vector3::new(10.0, 6.0, 4.0));
    /// # }
    /// ```
    #[inline]
    pub fn half_extents(&self) -> Vector<Real> {
        let half: Real = na::convert::<f64, Real>(0.5);
        (self.maxs - self.mins) * half
    }

    /// Returns the volume of this AABB.
    ///
    /// - **2D**: Returns the area (width × height)
    /// - **3D**: Returns the volume (width × height × depth)
    ///
    /// # Example
    ///
    /// ```rust
    /// # #[cfg(all(feature = "dim3", feature = "f32"))] {
    /// use parry3d::bounding_volume::Aabb;
    /// use nalgebra::Point3;
    ///
    /// // A 2x3x4 box
    /// let aabb = Aabb::new(
    ///     Point3::new(0.0, 0.0, 0.0),
    ///     Point3::new(2.0, 3.0, 4.0)
    /// );
    ///
    /// assert_eq!(aabb.volume(), 24.0); // 2 * 3 * 4
    /// # }
    /// ```
    #[inline]
    pub fn volume(&self) -> Real {
        let extents = self.extents();
        #[cfg(feature = "dim2")]
        return extents.x * extents.y;
        #[cfg(feature = "dim3")]
        return extents.x * extents.y * extents.z;
    }

    /// In 3D, returns the half-area. In 2D returns the half-perimeter of the AABB.
    pub fn half_area_or_perimeter(&self) -> Real {
        #[cfg(feature = "dim2")]
        return self.half_perimeter();
        #[cfg(feature = "dim3")]
        return self.half_area();
    }

    /// The half perimeter of this `Aabb`.
    #[cfg(feature = "dim2")]
    pub fn half_perimeter(&self) -> Real {
        let extents = self.extents();
        extents.x + extents.y
    }

    /// The half area of this `Aabb`.
    #[cfg(feature = "dim3")]
    pub fn half_area(&self) -> Real {
        let extents = self.extents();
        extents.x * (extents.y + extents.z) + extents.y * extents.z
    }

    /// The extents of this `Aabb`.
    #[inline]
    pub fn extents(&self) -> Vector<Real> {
        self.maxs - self.mins
    }

    /// Enlarges this `Aabb` so it also contains the point `pt`.
    pub fn take_point(&mut self, pt: Point<Real>) {
        self.mins = self.mins.coords.inf(&pt.coords).into();
        self.maxs = self.maxs.coords.sup(&pt.coords).into();
    }

    /// Computes the `Aabb` bounding `self` transformed by `m`.
    #[inline]
    pub fn transform_by(&self, m: &Isometry<Real>) -> Self {
        let ls_center = self.center();
        let center = m * ls_center;
        let ws_half_extents = m.absolute_transform_vector(&self.half_extents());

        Aabb::new(center + (-ws_half_extents), center + ws_half_extents)
    }

    /// Computes the Aabb bounding `self` translated by `translation`.
    #[inline]
    pub fn translated(mut self, translation: &Vector<Real>) -> Self {
        self.mins += translation;
        self.maxs += translation;
        self
    }

    #[inline]
    pub fn scaled(self, scale: &Vector<Real>) -> Self {
        let a = self.mins.coords.component_mul(scale);
        let b = self.maxs.coords.component_mul(scale);
        Self {
            mins: a.inf(&b).into(),
            maxs: a.sup(&b).into(),
        }
    }

    /// Returns an AABB with the same center as `self` but with extents scaled by `scale`.
    ///
    /// # Parameters
    /// - `scale`: the scaling factor. It can be non-uniform and/or negative. The AABB being
    ///   symmetric wrt. its center, a negative scale value has the same effect as scaling
    ///   by its absolute value.
    #[inline]
    #[must_use]
    pub fn scaled_wrt_center(self, scale: &Vector<Real>) -> Self {
        let center = self.center();
        // Multiply the extents by the scale. Negative scaling might modify the half-extent
        // sign, so we take the absolute value. The AABB being symmetric that absolute value
        // is  valid.
        let half_extents = self.half_extents().component_mul(scale).abs();
        Self::from_half_extents(center, half_extents)
    }

    /// The smallest bounding sphere containing this `Aabb`.
    #[inline]
    pub fn bounding_sphere(&self) -> BoundingSphere {
        let center = self.center();
        let radius = na::distance(&self.mins, &self.maxs) * 0.5;
        BoundingSphere::new(center, radius)
    }

    /// Does this AABB contains a point expressed in the same coordinate frame as `self`?
    #[inline]
    pub fn contains_local_point(&self, point: &Point<Real>) -> bool {
        for i in 0..DIM {
            if point[i] < self.mins[i] || point[i] > self.maxs[i] {
                return false;
            }
        }

        true
    }

    /// Computes the distance between the origin and this AABB.
    pub fn distance_to_origin(&self) -> Real {
        self.mins
            .coords
            .sup(&-self.maxs.coords)
            .sup(&Vector::zeros())
            .norm()
    }

    /// Does this AABB intersects an AABB `aabb2` moving at velocity `vel12` relative to `self`?
    #[inline]
    pub fn intersects_moving_aabb(&self, aabb2: &Self, vel12: Vector<Real>) -> bool {
        // Minkowski sum.
        let msum = Aabb {
            mins: self.mins - aabb2.maxs.coords,
            maxs: self.maxs - aabb2.mins.coords,
        };
        let ray = Ray::new(Point::origin(), vel12);

        msum.intersects_local_ray(&ray, 1.0)
    }

    /// Computes the intersection of this `Aabb` and another one.
    pub fn intersection(&self, other: &Aabb) -> Option<Aabb> {
        let result = Aabb {
            mins: Point::from(self.mins.coords.sup(&other.mins.coords)),
            maxs: Point::from(self.maxs.coords.inf(&other.maxs.coords)),
        };

        for i in 0..DIM {
            if result.mins[i] > result.maxs[i] {
                return None;
            }
        }

        Some(result)
    }

    /// Computes two AABBs for the intersection between two translated and rotated AABBs.
    ///
    /// This method returns two AABBs: the first is expressed in the local-space of `self`,
    /// and the second is expressed in the local-space of `aabb2`.
    pub fn aligned_intersections(
        &self,
        pos12: &Isometry<Real>,
        aabb2: &Self,
    ) -> Option<(Aabb, Aabb)> {
        let pos21 = pos12.inverse();

        let aabb2_1 = aabb2.transform_by(pos12);
        let inter1_1 = self.intersection(&aabb2_1)?;
        let inter1_2 = inter1_1.transform_by(&pos21);

        let aabb1_2 = self.transform_by(&pos21);
        let inter2_2 = aabb2.intersection(&aabb1_2)?;
        let inter2_1 = inter2_2.transform_by(pos12);

        Some((
            inter1_1.intersection(&inter2_1)?,
            inter1_2.intersection(&inter2_2)?,
        ))
    }

    /// Returns the difference between this `Aabb` and `rhs`.
    ///
    /// Removing another `Aabb` from `self` will result in zero, one, or up to 4 (in 2D) or 8 (in 3D)
    /// new smaller Aabbs.
    pub fn difference(&self, rhs: &Aabb) -> ArrayVec<Self, TWO_DIM> {
        self.difference_with_cut_sequence(rhs).0
    }

    /// Returns the difference between this `Aabb` and `rhs`.
    ///
    /// Removing another `Aabb` from `self` will result in zero, one, or up to 4 (in 2D) or 8 (in 3D)
    /// new smaller Aabbs.
    ///
    /// # Return
    /// This returns a pair where the first item are the new Aabbs and the second item is
    /// the sequence of cuts applied to `self` to obtain the new Aabbs. Each cut is performed
    /// along one axis identified by `-1, -2, -3` for `-X, -Y, -Z` and `1, 2, 3` for `+X, +Y, +Z`, and
    /// the plane’s bias.
    ///
    /// The cuts are applied sequentially. For example, if `result.1[0]` contains `1`, then it means
    /// that `result.0[0]` is equal to the piece of `self` lying in the negative half-space delimited
    /// by the plane with outward normal `+X`. Then, the other piece of `self` generated by this cut
    /// (i.e. the piece of `self` lying in the positive half-space delimited by the plane with outward
    /// normal `+X`) is the one that will be affected by the next cut.
    ///
    /// The returned cut sequence will be empty if the aabbs are disjoint.
    pub fn difference_with_cut_sequence(
        &self,
        rhs: &Aabb,
    ) -> (ArrayVec<Self, TWO_DIM>, ArrayVec<(i8, Real), TWO_DIM>) {
        let mut result = ArrayVec::new();
        let mut cut_sequence = ArrayVec::new();

        // NOTE: special case when the boxes are disjoint.
        //       This isn’t exactly the same as `!self.intersects(rhs)`
        //       because of the equality.
        for i in 0..DIM {
            if self.mins[i] >= rhs.maxs[i] || self.maxs[i] <= rhs.mins[i] {
                result.push(*self);
                return (result, cut_sequence);
            }
        }

        let mut rest = *self;

        for i in 0..DIM {
            if rhs.mins[i] > rest.mins[i] {
                let mut fragment = rest;
                fragment.maxs[i] = rhs.mins[i];
                rest.mins[i] = rhs.mins[i];
                result.push(fragment);
                cut_sequence.push((i as i8 + 1, rhs.mins[i]));
            }

            if rhs.maxs[i] < rest.maxs[i] {
                let mut fragment = rest;
                fragment.mins[i] = rhs.maxs[i];
                rest.maxs[i] = rhs.maxs[i];
                result.push(fragment);
                cut_sequence.push((-(i as i8 + 1), -rhs.maxs[i]));
            }
        }

        (result, cut_sequence)
    }

    /// Computes the vertices of this `Aabb`.
    ///
    /// The vertices are given in the following order in a right-handed coordinate system:
    /// ```text
    ///    y             3 - 2
    ///    |             |   |
    ///    ___ x         0 - 1
    /// ```
    #[inline]
    #[cfg(feature = "dim2")]
    pub fn vertices(&self) -> [Point<Real>; 4] {
        [
            Point::new(self.mins.x, self.mins.y),
            Point::new(self.maxs.x, self.mins.y),
            Point::new(self.maxs.x, self.maxs.y),
            Point::new(self.mins.x, self.maxs.y),
        ]
    }

    /// Computes the vertices of this `Aabb`.
    ///
    /// The vertices are given in the following order, in a right-handed coordinate system:
    /// ```text
    ///    y             3 - 2
    ///    |           7 − 6 |
    ///    ___ x       |   | 1  (the zero is below 3 and on the left of 1,
    ///   /            4 - 5     hidden by the 4-5-6-7 face.)
    ///  z
    /// ```
    #[inline]
    #[cfg(feature = "dim3")]
    pub fn vertices(&self) -> [Point<Real>; 8] {
        [
            Point::new(self.mins.x, self.mins.y, self.mins.z),
            Point::new(self.maxs.x, self.mins.y, self.mins.z),
            Point::new(self.maxs.x, self.maxs.y, self.mins.z),
            Point::new(self.mins.x, self.maxs.y, self.mins.z),
            Point::new(self.mins.x, self.mins.y, self.maxs.z),
            Point::new(self.maxs.x, self.mins.y, self.maxs.z),
            Point::new(self.maxs.x, self.maxs.y, self.maxs.z),
            Point::new(self.mins.x, self.maxs.y, self.maxs.z),
        ]
    }

    /// Splits this `Aabb` at its center, into four parts (as in a quad-tree).
    #[inline]
    #[cfg(feature = "dim2")]
    pub fn split_at_center(&self) -> [Aabb; 4] {
        let center = self.center();

        [
            Aabb::new(self.mins, center),
            Aabb::new(
                Point::new(center.x, self.mins.y),
                Point::new(self.maxs.x, center.y),
            ),
            Aabb::new(center, self.maxs),
            Aabb::new(
                Point::new(self.mins.x, center.y),
                Point::new(center.x, self.maxs.y),
            ),
        ]
    }

    /// Splits this `Aabb` at its center, into eight parts (as in an octree).
    #[inline]
    #[cfg(feature = "dim3")]
    pub fn split_at_center(&self) -> [Aabb; 8] {
        let center = self.center();

        [
            Aabb::new(
                Point::new(self.mins.x, self.mins.y, self.mins.z),
                Point::new(center.x, center.y, center.z),
            ),
            Aabb::new(
                Point::new(center.x, self.mins.y, self.mins.z),
                Point::new(self.maxs.x, center.y, center.z),
            ),
            Aabb::new(
                Point::new(center.x, center.y, self.mins.z),
                Point::new(self.maxs.x, self.maxs.y, center.z),
            ),
            Aabb::new(
                Point::new(self.mins.x, center.y, self.mins.z),
                Point::new(center.x, self.maxs.y, center.z),
            ),
            Aabb::new(
                Point::new(self.mins.x, self.mins.y, center.z),
                Point::new(center.x, center.y, self.maxs.z),
            ),
            Aabb::new(
                Point::new(center.x, self.mins.y, center.z),
                Point::new(self.maxs.x, center.y, self.maxs.z),
            ),
            Aabb::new(
                Point::new(center.x, center.y, center.z),
                Point::new(self.maxs.x, self.maxs.y, self.maxs.z),
            ),
            Aabb::new(
                Point::new(self.mins.x, center.y, center.z),
                Point::new(center.x, self.maxs.y, self.maxs.z),
            ),
        ]
    }

    /// Enlarges this AABB on each side by the given `half_extents`.
    #[must_use]
    pub fn add_half_extents(&self, half_extents: &Vector<Real>) -> Self {
        Self {
            mins: self.mins - half_extents,
            maxs: self.maxs + half_extents,
        }
    }

    /// Projects every point of `Aabb` on an arbitrary axis.
    pub fn project_on_axis(&self, axis: &UnitVector<Real>) -> (Real, Real) {
        let cuboid = Cuboid::new(self.half_extents());
        let shift = cuboid
            .local_support_point_toward(axis)
            .coords
            .dot(axis)
            .abs();
        let center = self.center().coords.dot(axis);
        (center - shift, center + shift)
    }

    #[cfg(feature = "dim3")]
    #[cfg(feature = "alloc")]
    pub fn intersects_spiral(
        &self,
        point: &Point<Real>,
        center: &Point<Real>,
        axis: &UnitVector<Real>,
        linvel: &Vector<Real>,
        angvel: Real,
    ) -> bool {
        use crate::utils::WBasis;
        use crate::utils::{Interval, IntervalFunction};
        use alloc::vec;

        struct SpiralPlaneDistance {
            center: Point<Real>,
            tangents: [Vector<Real>; 2],
            linvel: Vector<Real>,
            angvel: Real,
            point: na::Vector2<Real>,
            plane: Vector<Real>,
            bias: Real,
        }

        impl SpiralPlaneDistance {
            fn spiral_pt_at(&self, t: Real) -> Point<Real> {
                let angle = t * self.angvel;

                // NOTE: we construct the rotation matrix explicitly here instead
                //       of using `Rotation2::new()` because we will use similar
                //       formulas on the interval methods.
                let (sin, cos) = angle.sin_cos();
                let rotmat = na::Matrix2::new(cos, -sin, sin, cos);

                let rotated_pt = rotmat * self.point;
                let shift = self.tangents[0] * rotated_pt.x + self.tangents[1] * rotated_pt.y;
                self.center + self.linvel * t + shift
            }
        }

        impl IntervalFunction<Real> for SpiralPlaneDistance {
            fn eval(&self, t: Real) -> Real {
                let point_pos = self.spiral_pt_at(t);
                point_pos.coords.dot(&self.plane) - self.bias
            }

            fn eval_interval(&self, t: Interval<Real>) -> Interval<Real> {
                // This is the same as `self.eval` except that `t` is an interval.
                let angle = t * self.angvel;
                let (sin, cos) = angle.sin_cos();
                let rotmat = na::Matrix2::new(cos, -sin, sin, cos);

                let rotated_pt = rotmat * self.point.map(Interval::splat);
                let shift = self.tangents[0].map(Interval::splat) * rotated_pt.x
                    + self.tangents[1].map(Interval::splat) * rotated_pt.y;
                let point_pos =
                    self.center.map(Interval::splat) + self.linvel.map(Interval::splat) * t + shift;
                point_pos.coords.dot(&self.plane.map(Interval::splat)) - Interval::splat(self.bias)
            }

            fn eval_interval_gradient(&self, t: Interval<Real>) -> Interval<Real> {
                let angle = t * self.angvel;
                let (sin, cos) = angle.sin_cos();
                let rotmat = na::Matrix2::new(-sin, -cos, cos, -sin) * Interval::splat(self.angvel);

                let rotated_pt = rotmat * self.point.map(Interval::splat);
                let shift = self.tangents[0].map(Interval::splat) * rotated_pt.x
                    + self.tangents[1].map(Interval::splat) * rotated_pt.y;
                let point_vel = shift + self.linvel.map(Interval::splat);
                point_vel.dot(&self.plane.map(Interval::splat))
            }
        }

        let tangents = axis.orthonormal_basis();
        let dpos = point - center;
        let mut distance_fn = SpiralPlaneDistance {
            center: *center,
            tangents,
            linvel: *linvel,
            angvel,
            point: na::Vector2::new(dpos.dot(&tangents[0]), dpos.dot(&tangents[1])),
            plane: Vector::x(),
            bias: 0.0,
        };

        // Check the 8 planar faces of the Aabb.
        let mut roots = vec![];
        let mut candidates = vec![];

        let planes = [
            (-self.mins[0], -Vector::x(), 0),
            (self.maxs[0], Vector::x(), 0),
            (-self.mins[1], -Vector::y(), 1),
            (self.maxs[1], Vector::y(), 1),
            (-self.mins[2], -Vector::z(), 2),
            (self.maxs[2], Vector::z(), 2),
        ];

        let range = self.project_on_axis(axis);
        let range_bias = center.coords.dot(axis);
        let interval = Interval::sort(range.0, range.1) - range_bias;

        for (bias, axis, i) in &planes {
            distance_fn.plane = *axis;
            distance_fn.bias = *bias;

            crate::utils::find_root_intervals_to(
                &distance_fn,
                interval,
                1.0e-5,
                1.0e-5,
                100,
                &mut roots,
                &mut candidates,
            );

            for root in roots.drain(..) {
                let point = distance_fn.spiral_pt_at(root.midpoint());
                let (j, k) = ((i + 1) % 3, (i + 2) % 3);
                if point[j] >= self.mins[j]
                    && point[j] <= self.maxs[j]
                    && point[k] >= self.mins[k]
                    && point[k] <= self.maxs[k]
                {
                    return true;
                }
            }
        }

        false
    }
}

impl BoundingVolume for Aabb {
    #[inline]
    fn center(&self) -> Point<Real> {
        self.center()
    }

    #[inline]
    fn intersects(&self, other: &Aabb) -> bool {
        na::partial_le(&self.mins, &other.maxs) && na::partial_ge(&self.maxs, &other.mins)
    }

    #[inline]
    fn contains(&self, other: &Aabb) -> bool {
        na::partial_le(&self.mins, &other.mins) && na::partial_ge(&self.maxs, &other.maxs)
    }

    #[inline]
    fn merge(&mut self, other: &Aabb) {
        self.mins = self.mins.inf(&other.mins);
        self.maxs = self.maxs.sup(&other.maxs);
    }

    #[inline]
    fn merged(&self, other: &Aabb) -> Aabb {
        Aabb {
            mins: self.mins.inf(&other.mins),
            maxs: self.maxs.sup(&other.maxs),
        }
    }

    #[inline]
    fn loosen(&mut self, amount: Real) {
        assert!(amount >= 0.0, "The loosening margin must be positive.");
        self.mins += Vector::repeat(-amount);
        self.maxs += Vector::repeat(amount);
    }

    #[inline]
    fn loosened(&self, amount: Real) -> Aabb {
        assert!(amount >= 0.0, "The loosening margin must be positive.");
        Aabb {
            mins: self.mins + Vector::repeat(-amount),
            maxs: self.maxs + Vector::repeat(amount),
        }
    }

    #[inline]
    fn tighten(&mut self, amount: Real) {
        assert!(amount >= 0.0, "The tightening margin must be positive.");
        self.mins += Vector::repeat(amount);
        self.maxs += Vector::repeat(-amount);
        assert!(
            na::partial_le(&self.mins, &self.maxs),
            "The tightening margin is to large."
        );
    }

    #[inline]
    fn tightened(&self, amount: Real) -> Aabb {
        assert!(amount >= 0.0, "The tightening margin must be positive.");

        Aabb::new(
            self.mins + Vector::repeat(amount),
            self.maxs + Vector::repeat(-amount),
        )
    }
}
