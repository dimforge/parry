use crate::bounding_volume::AABB;
use crate::math::{Point, Real, SimdBool, SimdReal, Vector, DIM, SIMD_WIDTH};
use crate::query::SimdRay;
use crate::utils;
use num::{One, Zero};
use simba::simd::{SimdPartialOrd, SimdValue};

/// Four AABB represented as a single SoA AABB with SIMD components.
#[derive(Debug, Copy, Clone)]
pub struct SimdAABB {
    /// The min coordinates of the AABBs.
    pub mins: Point<SimdReal>,
    /// The max coordinates the AABBs.
    pub maxs: Point<SimdReal>,
}

#[cfg(feature = "serde-serialize")]
impl serde::Serialize for SimdAABB {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        use serde::ser::SerializeStruct;

        let mins: Point<[Real; SIMD_WIDTH]> = Point::from(
            self.mins
                .coords
                .map(|e| array![|ii| e.extract(ii); SIMD_WIDTH]),
        );
        let maxs: Point<[Real; SIMD_WIDTH]> = Point::from(
            self.maxs
                .coords
                .map(|e| array![|ii| e.extract(ii); SIMD_WIDTH]),
        );

        let mut simd_aabb = serializer.serialize_struct("simd_aabb", 2)?;
        simd_aabb.serialize_field("mins", &mins)?;
        simd_aabb.serialize_field("maxs", &maxs)?;
        simd_aabb.end()
    }
}

#[cfg(feature = "serde-serialize")]
impl<'de> serde::Deserialize<'de> for SimdAABB {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        struct Visitor {}

        impl<'de> serde::de::Visitor<'de> for Visitor {
            type Value = SimdAABB;
            fn expecting(&self, formatter: &mut std::fmt::Formatter) -> std::fmt::Result {
                write!(
                    formatter,
                    "two arrays containing at least {} floats",
                    SIMD_WIDTH * DIM * 2
                )
            }

            fn visit_seq<A>(self, mut seq: A) -> Result<Self::Value, A::Error>
            where
                A: serde::de::SeqAccess<'de>,
            {
                let mins: Point<[Real; SIMD_WIDTH]> = seq
                    .next_element()?
                    .ok_or_else(|| serde::de::Error::invalid_length(0, &self))?;
                let maxs: Point<[Real; SIMD_WIDTH]> = seq
                    .next_element()?
                    .ok_or_else(|| serde::de::Error::invalid_length(1, &self))?;
                let mins = Point::from(mins.coords.map(|e| SimdReal::from(e)));
                let maxs = Point::from(maxs.coords.map(|e| SimdReal::from(e)));
                Ok(SimdAABB { mins, maxs })
            }
        }

        deserializer.deserialize_struct("SimdAABB", &["mins", "maxs"], Visitor {})
    }
}

impl SimdAABB {
    /// An invalid AABB.
    pub fn new_invalid() -> Self {
        Self::splat(AABB::new_invalid())
    }

    /// Builds an SIMD aabb composed of four identical aabbs.
    pub fn splat(aabb: AABB) -> Self {
        Self {
            mins: Point::splat(aabb.mins),
            maxs: Point::splat(aabb.maxs),
        }
    }

    /// The center of all the AABBs represented by `self``.
    pub fn center(&self) -> Point<SimdReal> {
        na::center(&self.mins, &self.maxs)
    }

    /// The half-extents of all the AABBs represented by `self``.
    pub fn half_extents(&self) -> Vector<SimdReal> {
        (self.maxs - self.mins) * SimdReal::splat(0.5)
    }

    /// The radius of all the AABBs represented by `self``.
    pub fn radius(&self) -> SimdReal {
        (self.maxs - self.mins).norm()
    }

    /// Dilate all the AABBs represented by `self`` by their extents multiplied
    /// by the given scale `factor`.
    pub fn dilate_by_factor(&mut self, factor: SimdReal) {
        // If some of the AABBs on this SimdAABB are invalid,
        // don't, dilate them.
        let is_valid = self.mins.x.simd_le(self.maxs.x);
        let factor = factor.select(is_valid, SimdReal::zero());

        // NOTE: we multiply each by factor instead of doing
        // (maxs - mins) * factor. That's to avoid overflows (and
        // therefore NaNs if this SimdAABB contains some invalid
        // AABBs initialised with Real::MAX
        let dilation = self.maxs * factor - self.mins * factor;
        self.mins -= dilation;
        self.maxs += dilation;
    }

    /// Replace the `i-th` AABB of this SIMD AAAB by the given value.
    pub fn replace(&mut self, i: usize, aabb: AABB) {
        self.mins.replace(i, aabb.mins);
        self.maxs.replace(i, aabb.maxs);
    }

    /// Casts a ray on all the AABBs represented by `self`.
    pub fn cast_local_ray(&self, ray: &SimdRay, max_toi: SimdReal) -> (SimdBool, SimdReal) {
        let zero = SimdReal::zero();
        let one = SimdReal::one();
        let infinity = SimdReal::splat(Real::MAX);

        let mut hit = SimdBool::splat(true);
        let mut tmin = SimdReal::zero();
        let mut tmax = max_toi;

        // TODO: could this be optimized more considering we really just need a boolean answer?
        for i in 0usize..DIM {
            let is_not_zero = ray.dir[i].simd_ne(zero);
            let is_zero_test =
                ray.origin[i].simd_ge(self.mins[i]) & ray.origin[i].simd_le(self.maxs[i]);
            let is_not_zero_test = {
                let denom = one / ray.dir[i];
                let mut inter_with_near_plane =
                    ((self.mins[i] - ray.origin[i]) * denom).select(is_not_zero, -infinity);
                let mut inter_with_far_plane =
                    ((self.maxs[i] - ray.origin[i]) * denom).select(is_not_zero, infinity);

                let gt = inter_with_near_plane.simd_gt(inter_with_far_plane);
                utils::simd_swap(gt, &mut inter_with_near_plane, &mut inter_with_far_plane);

                tmin = tmin.simd_max(inter_with_near_plane);
                tmax = tmax.simd_min(inter_with_far_plane);

                tmin.simd_le(tmax)
            };

            hit = hit & is_not_zero_test.select(is_not_zero, is_zero_test);
        }

        (hit, tmin)
    }

    /// Computes the distances between a point and all the AABBs represented by `self`.
    pub fn distance_to_local_point(&self, point: &Point<SimdReal>) -> SimdReal {
        let mins_point = self.mins - point;
        let point_maxs = point - self.maxs;
        let shift = mins_point.sup(&point_maxs).sup(&na::zero());
        shift.norm()
    }

    /// Computes the distances between the origin and all the AABBs represented by `self`.
    pub fn distance_to_origin(&self) -> SimdReal {
        self.mins
            .coords
            .sup(&-self.maxs.coords)
            .sup(&Vector::zeros())
            .norm()
    }

    /// Check which AABB represented by `self` contains the given `point`.
    pub fn contains_local_point(&self, point: &Point<SimdReal>) -> SimdBool {
        #[cfg(feature = "dim2")]
        return self.mins.x.simd_le(point.x)
            & self.mins.y.simd_le(point.y)
            & self.maxs.x.simd_ge(point.x)
            & self.maxs.y.simd_ge(point.y);

        #[cfg(feature = "dim3")]
        return self.mins.x.simd_le(point.x)
            & self.mins.y.simd_le(point.y)
            & self.mins.z.simd_le(point.z)
            & self.maxs.x.simd_ge(point.x)
            & self.maxs.y.simd_ge(point.y)
            & self.maxs.z.simd_ge(point.z);
    }

    /// Lanewise check which AABB represented by `self` contains the given set of `other` aabbs.
    /// The check is performed lane-wise.
    #[cfg(feature = "dim2")]
    pub fn contains(&self, other: &SimdAABB) -> SimdBool {
        self.mins.x.simd_le(other.mins.x)
            & self.mins.y.simd_le(other.mins.y)
            & self.maxs.x.simd_ge(other.maxs.x)
            & self.maxs.y.simd_ge(other.maxs.y)
    }

    /// Lanewise check which AABB represented by `self` contains the given set of `other` aabbs.
    /// The check is performed lane-wise.
    #[cfg(feature = "dim3")]
    pub fn contains(&self, other: &SimdAABB) -> SimdBool {
        self.mins.x.simd_le(other.mins.x)
            & self.mins.y.simd_le(other.mins.y)
            & self.mins.z.simd_le(other.mins.z)
            & self.maxs.x.simd_ge(other.maxs.x)
            & self.maxs.y.simd_ge(other.maxs.y)
            & self.maxs.z.simd_ge(other.maxs.z)
    }

    /// Lanewise check which AABB represented by `self` intersects the given set of `other` aabbs.
    /// The check is performed lane-wise.
    #[cfg(feature = "dim2")]
    pub fn intersects(&self, other: &SimdAABB) -> SimdBool {
        self.mins.x.simd_le(other.maxs.x)
            & other.mins.x.simd_le(self.maxs.x)
            & self.mins.y.simd_le(other.maxs.y)
            & other.mins.y.simd_le(self.maxs.y)
    }

    /// Check which AABB represented by `self` contains the given set of `other` aabbs.
    /// The check is performed lane-wise.
    #[cfg(feature = "dim3")]
    pub fn intersects(&self, other: &SimdAABB) -> SimdBool {
        self.mins.x.simd_le(other.maxs.x)
            & other.mins.x.simd_le(self.maxs.x)
            & self.mins.y.simd_le(other.maxs.y)
            & other.mins.y.simd_le(self.maxs.y)
            & self.mins.z.simd_le(other.maxs.z)
            & other.mins.z.simd_le(self.maxs.z)
    }

    /// Merge all the AABB represented by `self` into a single one.
    pub fn to_merged_aabb(&self) -> AABB {
        AABB::new(
            self.mins.coords.map(|e| e.simd_horizontal_min()).into(),
            self.maxs.coords.map(|e| e.simd_horizontal_max()).into(),
        )
    }
}

impl From<[AABB; SIMD_WIDTH]> for SimdAABB {
    fn from(aabbs: [AABB; SIMD_WIDTH]) -> Self {
        let mins = array![|ii| aabbs[ii].mins; SIMD_WIDTH];
        let maxs = array![|ii| aabbs[ii].maxs; SIMD_WIDTH];

        SimdAABB {
            mins: Point::from(mins),
            maxs: Point::from(maxs),
        }
    }
}
