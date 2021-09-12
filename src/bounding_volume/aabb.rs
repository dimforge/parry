//! Axis Aligned Bounding Box.

use crate::bounding_volume::{BoundingSphere, BoundingVolume};
use crate::math::{Isometry, Point, Real, UnitVector, Vector, DIM};
use crate::shape::{Cuboid, SupportMap};
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
    /// The vertex indices of each edge of this AABB.
    ///
    /// This gives, for each edge of this AABB, the indices of its
    /// vertices when taken from the `self.vertices()` array.
    /// Here is how the faces are numbered, assuming
    /// a right-handed coordinate system:
    ///
    ///    y             3 - 2
    ///    |           7 − 6 |
    ///    ___ x       |   | 1  (the zero is bellow 3 and on the left of 1, hidden by the 4-5-6-7 face.)
    ///   /            4 - 5
    ///  z
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

    /// The vertex indices of each face of this AABB.
    ///
    /// This gives, for each face of this AABB, the indices of its
    /// vertices when taken from the `self.vertices()` array.
    /// Here is how the faces are numbered, assuming
    /// a right-handed coordinate system:
    ///
    ///    y             3 - 2
    ///    |           7 − 6 |
    ///    ___ x       |   | 1  (the zero is bellow 3 and on the left of 1, hidden by the 4-5-6-7 face.)
    ///   /            4 - 5
    ///  z
    #[cfg(feature = "dim3")]
    pub const FACES_VERTEX_IDS: [(usize, usize, usize, usize); 6] = [
        (1, 2, 6, 5),
        (0, 3, 7, 4),
        (2, 3, 7, 6),
        (1, 0, 4, 5),
        (4, 5, 6, 7),
        (0, 1, 2, 3),
    ];

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

    /// The volume of this AABB.
    #[inline]
    pub fn volume(&self) -> Real {
        let extents = self.extents();
        #[cfg(feature = "dim2")]
        return extents.x * extents.y;
        #[cfg(feature = "dim3")]
        return extents.x * extents.y * extents.z;
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
            Point::new(self.maxs.x, self.mins.y, self.mins.z),
            Point::new(self.maxs.x, self.maxs.y, self.mins.z),
            Point::new(self.mins.x, self.maxs.y, self.mins.z),
            Point::new(self.mins.x, self.mins.y, self.maxs.z),
            Point::new(self.maxs.x, self.mins.y, self.maxs.z),
            Point::new(self.maxs.x, self.maxs.y, self.maxs.z),
            Point::new(self.mins.x, self.maxs.y, self.maxs.z),
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

    /// Projects every point of AABB on an arbitrary axis.
    pub fn project_on_axis(&self, axis: &UnitVector<Real>) -> (Real, Real) {
        let cuboid = Cuboid::new(self.half_extents());
        let shift = cuboid
            .local_support_point_toward(axis)
            .coords
            .dot(&axis)
            .abs();
        let center = self.center().coords.dot(&axis);
        (center - shift, center + shift)
    }

    #[cfg(feature = "dim3")]
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
                //       formulaes on the interval methods.
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

        // Check the 8 planar faces of the AABB.
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

        let range = self.project_on_axis(&axis);
        let range_bias = center.coords.dot(&axis);
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
