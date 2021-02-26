//! Axis Aligned Bounding Box.

use crate::bounding_volume::{BoundingSphere, BoundingVolume};
use crate::math::{Isometry, Point, Real, Vector, DIM};
use crate::utils::IsometryOps;
use na;
use num::Bounded;

/// An Axis Aligned Bounding Box.
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[derive(Debug, PartialEq, Copy, Clone)]
pub struct AABB {
    pub mins: Point<Real>,
    pub maxs: Point<Real>,
}

impl AABB {
    /// Creates a new AABB.
    ///
    /// # Arguments:
    ///   * `mins` - position of the point with the smallest coordinates.
    ///   * `maxs` - position of the point with the highest coordinates. Each component of `mins`
    ///   must be smaller than the related components of `maxs`.
    #[inline]
    pub fn new(mins: Point<Real>, maxs: Point<Real>) -> AABB {
        AABB { mins, maxs }
    }

    /// Creates an invalid AABB with `mins` components set to `Real::max_values` and `maxs`components set to `-Real::max_values`.
    ///
    /// This is often used as the initial values of some AABB merging algorithms.
    #[inline]
    pub fn new_invalid() -> Self {
        Self::new(
            Vector::repeat(Real::max_value()).into(),
            Vector::repeat(-Real::max_value()).into(),
        )
    }

    /// Creates a new AABB from its center and its half-extents.
    #[inline]
    pub fn from_half_extents(center: Point<Real>, half_extents: Vector<Real>) -> Self {
        Self::new(center - half_extents, center + half_extents)
    }

    /// Creates a new AABB from a set of points.
    pub fn from_points<'a, I>(pts: I) -> Self
    where
        I: IntoIterator<Item = &'a Point<Real>>,
    {
        super::aabb_utils::local_point_cloud_aabb(pts)
    }

    /// The center of this AABB.
    #[inline]
    pub fn center(&self) -> Point<Real> {
        na::center(&self.mins, &self.maxs)
    }

    /// The half extents of this AABB.
    #[inline]
    pub fn half_extents(&self) -> Vector<Real> {
        let half: Real = na::convert::<f64, Real>(0.5);
        (self.maxs - self.mins) * half
    }

    /// The extents of this AABB.
    #[inline]
    pub fn extents(&self) -> Vector<Real> {
        self.maxs - self.mins
    }

    /// Enlarges this AABB so it also contains the point `pt`.
    pub fn take_point(&mut self, pt: Point<Real>) {
        self.mins = self.mins.coords.inf(&pt.coords).into();
        self.maxs = self.maxs.coords.sup(&pt.coords).into();
    }

    /// Computes the AABB bounding `self` transformed by `m`.
    #[inline]
    pub fn transform_by(&self, m: &Isometry<Real>) -> Self {
        let ls_center = self.center();
        let center = m * ls_center;
        let ws_half_extents = m.absolute_transform_vector(&self.half_extents());

        AABB::new(center + (-ws_half_extents), center + ws_half_extents)
    }

    /// The smallest bounding sphere containing this AABB.
    #[inline]
    pub fn bounding_sphere(&self) -> BoundingSphere {
        let center = self.center();
        let rad = na::distance(&self.mins, &self.maxs);

        BoundingSphere::new(center, rad)
    }

    #[inline]
    pub fn contains_local_point(&self, point: &Point<Real>) -> bool {
        for i in 0..DIM {
            if point[i] < self.mins[i] || point[i] > self.maxs[i] {
                return false;
            }
        }

        true
    }

    /// Computes the vertices of this AABB.
    #[inline]
    #[cfg(feature = "dim2")]
    pub fn vertices(&self) -> [Point<Real>; 4] {
        [
            Point::new(self.mins.x, self.mins.y),
            Point::new(self.mins.x, self.maxs.y),
            Point::new(self.maxs.x, self.mins.y),
            Point::new(self.maxs.x, self.maxs.y),
        ]
    }

    /// Computes the vertices of this AABB.
    #[inline]
    #[cfg(feature = "dim3")]
    pub fn vertices(&self) -> [Point<Real>; 8] {
        [
            Point::new(self.mins.x, self.mins.y, self.mins.z),
            Point::new(self.mins.x, self.mins.y, self.maxs.z),
            Point::new(self.mins.x, self.maxs.y, self.mins.z),
            Point::new(self.mins.x, self.maxs.y, self.maxs.z),
            Point::new(self.maxs.x, self.mins.y, self.mins.z),
            Point::new(self.maxs.x, self.mins.y, self.maxs.z),
            Point::new(self.maxs.x, self.maxs.y, self.mins.z),
            Point::new(self.maxs.x, self.maxs.y, self.maxs.z),
        ]
    }

    /// Splits this AABB at its center, into four parts (as in a quad-tree).
    #[inline]
    #[cfg(feature = "dim2")]
    pub fn split_at_center(&self) -> [AABB; 4] {
        let center = self.center();

        [
            AABB::new(self.mins, center),
            AABB::new(
                Point::new(center.x, self.mins.y),
                Point::new(self.maxs.x, center.y),
            ),
            AABB::new(center, self.maxs),
            AABB::new(
                Point::new(self.mins.x, center.y),
                Point::new(center.x, self.maxs.y),
            ),
        ]
    }

    /// Splits this AABB at its center, into height parts (as in an octree).
    #[inline]
    #[cfg(feature = "dim3")]
    pub fn split_at_center(&self) -> [AABB; 8] {
        let center = self.center();

        [
            AABB::new(
                Point::new(self.mins.x, self.mins.y, self.mins.z),
                Point::new(center.x, center.y, center.z),
            ),
            AABB::new(
                Point::new(center.x, self.mins.y, self.mins.z),
                Point::new(self.maxs.x, center.y, center.z),
            ),
            AABB::new(
                Point::new(center.x, center.y, self.mins.z),
                Point::new(self.maxs.x, self.maxs.y, center.z),
            ),
            AABB::new(
                Point::new(self.mins.x, center.y, self.mins.z),
                Point::new(center.x, self.maxs.y, center.z),
            ),
            AABB::new(
                Point::new(self.mins.x, self.mins.y, center.z),
                Point::new(center.x, center.y, self.maxs.z),
            ),
            AABB::new(
                Point::new(center.x, self.mins.y, center.z),
                Point::new(self.maxs.x, center.y, self.maxs.z),
            ),
            AABB::new(
                Point::new(center.x, center.y, center.z),
                Point::new(self.maxs.x, self.maxs.y, self.maxs.z),
            ),
            AABB::new(
                Point::new(self.mins.x, center.y, center.z),
                Point::new(center.x, self.maxs.y, self.maxs.z),
            ),
        ]
    }
}

impl BoundingVolume for AABB {
    #[inline]
    fn center(&self) -> Point<Real> {
        self.center()
    }

    #[inline]
    fn intersects(&self, other: &AABB) -> bool {
        na::partial_le(&self.mins, &other.maxs) && na::partial_ge(&self.maxs, &other.mins)
    }

    #[inline]
    fn contains(&self, other: &AABB) -> bool {
        na::partial_le(&self.mins, &other.mins) && na::partial_ge(&self.maxs, &other.maxs)
    }

    #[inline]
    fn merge(&mut self, other: &AABB) {
        self.mins = self.mins.inf(&other.mins);
        self.maxs = self.maxs.sup(&other.maxs);
    }

    #[inline]
    fn merged(&self, other: &AABB) -> AABB {
        AABB {
            mins: self.mins.inf(&other.mins),
            maxs: self.maxs.sup(&other.maxs),
        }
    }

    #[inline]
    fn loosen(&mut self, amount: Real) {
        assert!(amount >= 0.0, "The loosening margin must be positive.");
        self.mins = self.mins + Vector::repeat(-amount);
        self.maxs = self.maxs + Vector::repeat(amount);
    }

    #[inline]
    fn loosened(&self, amount: Real) -> AABB {
        assert!(amount >= 0.0, "The loosening margin must be positive.");
        AABB {
            mins: self.mins + Vector::repeat(-amount),
            maxs: self.maxs + Vector::repeat(amount),
        }
    }

    #[inline]
    fn tighten(&mut self, amount: Real) {
        assert!(amount >= 0.0, "The tightening margin must be positive.");
        self.mins = self.mins + Vector::repeat(amount);
        self.maxs = self.maxs + Vector::repeat(-amount);
        assert!(
            na::partial_le(&self.mins, &self.maxs),
            "The tightening margin is to large."
        );
    }

    #[inline]
    fn tightened(&self, amount: Real) -> AABB {
        assert!(amount >= 0.0, "The tightening margin must be positive.");

        AABB::new(
            self.mins + Vector::repeat(amount),
            self.maxs + Vector::repeat(-amount),
        )
    }
}
