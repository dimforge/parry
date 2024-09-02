//! Traits and structure needed to cast rays.

use crate::math::{Isometry, Point, Real, Vector};
use crate::shape::FeatureId;

#[cfg(feature = "rkyv")]
use rkyv::{bytecheck, CheckBytes};

/// A Ray.
#[derive(Debug, Clone, Copy)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Deserialize, rkyv::Serialize, CheckBytes),
    archive(as = "Self")
)]
#[repr(C)]
pub struct Ray {
    /// Starting point of the ray.
    pub origin: Point<Real>,
    /// Direction of the ray.
    pub dir: Vector<Real>,
}

impl Ray {
    /// Creates a new ray starting from `origin` and with the direction `dir`.
    pub fn new(origin: Point<Real>, dir: Vector<Real>) -> Ray {
        Ray { origin, dir }
    }

    /// Transforms this ray by the given isometry.
    #[inline]
    pub fn transform_by(&self, m: &Isometry<Real>) -> Self {
        Self::new(m * self.origin, m * self.dir)
    }

    /// Transforms this ray by the inverse of the given isometry.
    #[inline]
    pub fn inverse_transform_by(&self, m: &Isometry<Real>) -> Self {
        Self::new(
            m.inverse_transform_point(&self.origin),
            m.inverse_transform_vector(&self.dir),
        )
    }

    /// Translates this ray by the given vector. Its direction is left unchanged.
    #[inline]
    pub fn translate_by(&self, v: Vector<Real>) -> Self {
        Self::new(self.origin + v, self.dir)
    }

    /// Computes the point at the given parameter on this line.
    ///
    /// This computes `self.origin + self.dir * t`.
    #[inline]
    pub fn point_at(&self, t: Real) -> Point<Real> {
        self.origin + self.dir * t
    }
}

/// Structure containing the result of a successful ray cast.
#[derive(Copy, Clone, Debug)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Deserialize, rkyv::Serialize, CheckBytes),
    archive(as = "Self")
)]
pub struct RayIntersection {
    /// The time of impact of the ray with the object. The exact contact point can be computed
    /// with: `ray.point_at(time_of_impact)` or equivalently `origin + dir * time_of_impact` where `origin` is the origin of the ray;
    /// `dir` is its direction and `time_of_impact` is the value of this field.
    pub time_of_impact: Real,

    /// The normal at the intersection point.
    ///
    /// If the origin of the ray is inside the shape and the shape is not solid,
    /// the normal will point towards the interior of the shape.
    /// Otherwise, the normal points outward.
    ///
    /// If the `time_of_impact` is exactly zero, the normal might not be reliable.
    // TODO: use a Unit<Vector> instead.
    pub normal: Vector<Real>,

    /// Feature at the intersection point.
    pub feature: FeatureId,
}

impl RayIntersection {
    #[inline]
    /// Creates a new `RayIntersection`.
    #[cfg(feature = "dim3")]
    pub fn new(time_of_impact: Real, normal: Vector<Real>, feature: FeatureId) -> RayIntersection {
        RayIntersection {
            time_of_impact,
            normal,
            feature,
        }
    }

    #[inline]
    /// Creates a new `RayIntersection`.
    #[cfg(feature = "dim2")]
    pub fn new(time_of_impact: Real, normal: Vector<Real>, feature: FeatureId) -> RayIntersection {
        RayIntersection {
            time_of_impact,
            normal,
            feature,
        }
    }

    #[inline]
    pub fn transform_by(&self, transform: &Isometry<Real>) -> Self {
        RayIntersection {
            time_of_impact: self.time_of_impact,
            normal: transform * self.normal,
            feature: self.feature,
        }
    }
}

/// Traits of objects which can be transformed and tested for intersection with a ray.
pub trait RayCast {
    /// Computes the time of impact between this transform shape and a ray.
    fn cast_local_ray(&self, ray: &Ray, max_time_of_impact: Real, solid: bool) -> Option<Real> {
        self.cast_local_ray_and_get_normal(ray, max_time_of_impact, solid)
            .map(|inter| inter.time_of_impact)
    }

    /// Computes the time of impact, and normal between this transformed shape and a ray.
    fn cast_local_ray_and_get_normal(
        &self,
        ray: &Ray,
        max_time_of_impact: Real,
        solid: bool,
    ) -> Option<RayIntersection>;

    /// Tests whether a ray intersects this transformed shape.
    #[inline]
    fn intersects_local_ray(&self, ray: &Ray, max_time_of_impact: Real) -> bool {
        self.cast_local_ray(ray, max_time_of_impact, true).is_some()
    }

    /// Computes the time of impact between this transform shape and a ray.
    fn cast_ray(
        &self,
        m: &Isometry<Real>,
        ray: &Ray,
        max_time_of_impact: Real,
        solid: bool,
    ) -> Option<Real> {
        let ls_ray = ray.inverse_transform_by(m);
        self.cast_local_ray(&ls_ray, max_time_of_impact, solid)
    }

    /// Computes the time of impact, and normal between this transformed shape and a ray.
    fn cast_ray_and_get_normal(
        &self,
        m: &Isometry<Real>,
        ray: &Ray,
        max_time_of_impact: Real,
        solid: bool,
    ) -> Option<RayIntersection> {
        let ls_ray = ray.inverse_transform_by(m);
        self.cast_local_ray_and_get_normal(&ls_ray, max_time_of_impact, solid)
            .map(|inter| inter.transform_by(m))
    }

    /// Tests whether a ray intersects this transformed shape.
    #[inline]
    fn intersects_ray(&self, m: &Isometry<Real>, ray: &Ray, max_time_of_impact: Real) -> bool {
        let ls_ray = ray.inverse_transform_by(m);
        self.intersects_local_ray(&ls_ray, max_time_of_impact)
    }
}
