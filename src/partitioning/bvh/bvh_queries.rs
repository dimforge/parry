use super::{Bvh, BvhNode};
use crate::bounding_volume::{Aabb, BoundingVolume};
use crate::math::Point;
use crate::math::Real;
use crate::query::PointProjection;
use crate::query::{PointQuery, Ray};
use crate::shape::FeatureId;

#[cfg(all(feature = "simd-is-enabled", feature = "dim3", feature = "f32"))]
pub(super) struct SimdInvRay {
    // TODO: we need to use `glam` here instead of `wide` because `wide` is lacking
    //       operations for getting the min/max vector element.
    //       We can switch back to `wide` once it's supported.
    pub origin: glam::Vec3A,
    pub inv_dir: glam::Vec3A,
}

#[cfg(all(feature = "simd-is-enabled", feature = "dim3", feature = "f32"))]
impl From<Ray> for SimdInvRay {
    fn from(ray: Ray) -> Self {
        let inv_dir = ray.dir.map(|r| {
            if r.abs() < Real::EPSILON {
                r.signum() / Real::EPSILON
            } else {
                1.0 / r
            }
        });
        Self {
            origin: glam::Vec3A::from([ray.origin.x, ray.origin.y, ray.origin.z]),
            inv_dir: glam::Vec3A::from([inv_dir.x, inv_dir.y, inv_dir.z]),
        }
    }
}

impl Bvh {
    /// Iterates through all the leaves with an AABB intersecting the given `aabb`.
    pub fn intersect_aabb<'a>(&'a self, aabb: &'a Aabb) -> impl Iterator<Item = u32> + 'a {
        self.leaves(|node: &BvhNode| node.aabb().intersects(aabb))
    }

    /// Projects a point on this BVH using the provided leaf projection function.
    ///
    /// The `primitive_check` delegates the point-projection task to an external function that
    /// is assumed to map a leaf index to an actual geometry to project on. The `Real` argument
    /// given to that closure is the distance to the closest point found so far (or is equal to
    /// `max_distance` if no projection was found so far).
    pub fn project_point(
        &self,
        point: &Point<Real>,
        max_distance: Real,
        primitive_check: impl Fn(u32, Real) -> Option<PointProjection>,
    ) -> Option<(u32, (Real, PointProjection))> {
        self.find_best(
            max_distance,
            |node: &BvhNode, _| node.aabb().distance_to_local_point(point, true),
            |primitive, _| {
                let proj = primitive_check(primitive, max_distance)?;
                Some((na::distance(&proj.point, point), proj))
            },
        )
    }

    /// Projects a point on this BVH using the provided leaf projection function.
    ///
    /// Also returns the feature the point was projected on.
    ///
    /// The `primitive_check` delegates the point-projection task to an external function that
    /// is assumed to map a leaf index to an actual geometry to project on. The `Real` argument
    /// given to that closure is the distance to the closest point found so far (or is equal to
    /// `max_distance` if no projection was found so far).
    pub fn project_point_and_get_feature(
        &self,
        point: &Point<Real>,
        max_distance: Real,
        primitive_check: impl Fn(u32, Real) -> Option<(PointProjection, FeatureId)>,
    ) -> Option<(u32, (Real, (PointProjection, FeatureId)))> {
        self.find_best(
            max_distance,
            |node: &BvhNode, _| node.aabb().distance_to_local_point(point, true),
            |primitive, _| {
                let proj = primitive_check(primitive, max_distance)?;
                Some((na::distance(&proj.0.point, point), proj))
            },
        )
    }

    /// Casts a ray on this BVH using the provided leaf ray-cast function.
    ///
    /// The `primitive_check` delegates the ray-casting task to an external function that
    /// is assumed to map a leaf index to an actual geometry to cast a ray on. The `Real` argument
    /// given to that closure is the distance to the closest ray hit found so far (or is equal to
    /// `max_time_of_impact` if no projection was found so far).
    #[cfg(not(all(feature = "simd-is-enabled", feature = "dim3", feature = "f32")))]
    pub fn cast_ray(
        &self,
        ray: &Ray,
        max_time_of_impact: Real,
        primitive_check: impl Fn(u32, Real) -> Option<Real>,
    ) -> Option<(u32, Real)> {
        self.find_best(
            max_time_of_impact,
            |node: &BvhNode, best_so_far| node.cast_ray(ray, best_so_far),
            primitive_check,
        )
    }

    /// Casts a ray on this BVH using the provided leaf ray-cast function.
    ///
    /// The `primitive_check` delegates the ray-casting task to an external function that
    /// is assumed to map a leaf index to an actual geometry to cast a ray on. The `Real` argument
    /// given to that closure is the distance to the closest ray hit found so far (or is equal to
    /// `max_time_of_impact` if no projection was found so far).
    #[cfg(all(feature = "simd-is-enabled", feature = "dim3", feature = "f32"))]
    pub fn cast_ray(
        &self,
        ray: &Ray,
        max_time_of_impact: Real,
        primitive_check: impl Fn(u32, Real) -> Option<Real>,
    ) -> Option<(u32, Real)> {
        // The commented code below relies on depth-first traversal instead of
        // depth first. Interesting to compare both approaches (the best-first is
        // expected to be consistently faster though).
        //
        // let simd_inv_ray = SimdInvRay::from(*ray);
        // let mut best_primitive = u32::MAX;
        // let mut best_toi = Real::MAX;
        // self.traverse(|node: &BvhNode| {
        //     if node.cast_inv_ray_simd(&simd_inv_ray) >= best_toi {
        //         return TraversalAction::Prune;
        //     }
        //
        //     if let Some(primitive) = node.leaf_data() {
        //         let toi = primitive_check(primitive as u32, best_toi);
        //         if toi < best_toi {
        //             best_primitive = primitive as u32;
        //             best_toi = toi;
        //         }
        //     }
        //
        //     TraversalAction::Continue
        // });
        // (best_primitive, best_toi)

        let simd_inv_ray = SimdInvRay::from(*ray);
        self.find_best(
            max_time_of_impact,
            |node: &BvhNode, _best_so_far| node.cast_inv_ray_simd(&simd_inv_ray),
            |primitive, best_so_far| primitive_check(primitive, best_so_far),
        )
    }
}
