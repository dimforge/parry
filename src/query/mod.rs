//! Non-persistent geometric queries.
//!
//! # General cases
//! The most general methods provided by this module are:
//!
//! * [`query::closest_points()`] to compute the closest points between two shapes.
//! * [`query::distance()`] to compute the distance between two shapes.
//! * [`query::contact()`] to compute one pair of contact points between two shapes, including penetrating contact.
//! * [`query::intersection_test()`] to determine if two shapes are intersecting or not.
//! * [`query::time_of_impact()`] to determine when two shapes undergoing translational motions hit for the first time.
//! * [`query::nonlinear_time_of_impact()`] to determine when two shapes undergoing continuous rigid motions hit for the first time.
//!
//! Ray-casting and point-projection can be achieved by importing traits:
//!
//! * [`query::RayCast`] for ray-casting.
//! * [`query::PointQuery`] for point projection.
//!
//! # Specific cases
//! The functions exported by the `details` submodule are more specific versions of the ones described above.
//! For example `distance_ball_ball` computes the distance between two shapes known at compile-time to be balls.
//! They are less convenient to use than the most generic version but will be slightly faster due to the lack of dynamic dispatch.
//! The specific functions have the form `[operation]_[shape1]_[shape2]()` where:
//!
//! * `[operation]` can be `closest_points`, `distance`, `contact`, `intersection_test` or `time_of_impact`.
//! * `[shape1]` is the type of the first shape passed to the function, e.g., `ball`, or `halfspace`. Can also identify a trait implemented by supported shapes, e.g., `support_map`.
//! * `[shape2]` is the type of the second shape passed to the function, e.g., `ball`, or `halfspace`. Can also identify a trait implemented by supported shapes, e.g., `support_map`.

pub use self::closest_points::{closest_points, ClosestPoints};
pub use self::contact::{contact, Contact};
pub use self::contact_manifolds::{
    ContactManifold, ContactManifoldsWorkspace, KinematicsCategory, TrackedContact,
};
pub use self::default_query_dispatcher::DefaultQueryDispatcher;
pub use self::distance::distance;
pub use self::error::Unsupported;
pub use self::intersection_test::intersection_test;
pub use self::nonlinear_time_of_impact::nonlinear_time_of_impact;
pub use self::point::{PointProjection, PointQuery, PointQueryWithLocation};
pub use self::query_dispatcher::{
    PersistentQueryDispatcher, QueryDispatcher, QueryDispatcherChain,
};
pub use self::ray::{Ray, RayCast, RayIntersection, SimdRay};
pub use self::time_of_impact::{time_of_impact, TOIStatus, TOI};

mod clip;
pub mod closest_points;
pub mod contact;
mod contact_manifolds;
mod default_query_dispatcher;
mod distance;
pub mod epa;
mod error;
pub mod gjk;
mod intersection_test;
mod nonlinear_time_of_impact;
pub mod point;
mod query_dispatcher;
mod ray;
pub mod sat;
mod time_of_impact;
pub mod visitors;

/// Queries dedicated to specific pairs of shapes.
pub mod details {
    pub use super::clip::*;
    pub use super::closest_points::*;
    pub use super::contact::{
        contact_ball_ball, contact_ball_convex_polyhedron, contact_composite_shape_shape,
        contact_convex_polyhedron_ball, contact_halfspace_support_map,
        contact_shape_composite_shape, contact_support_map_halfspace,
        contact_support_map_support_map, contact_support_map_support_map_with_params,
    };
    pub use super::distance::{
        distance_ball_ball, distance_composite_shape_shape, distance_halfspace_support_map,
        distance_shape_composite_shape, distance_support_map_halfspace,
        distance_support_map_support_map, distance_support_map_support_map_with_params,
    };
    pub use super::intersection_test::*;
    pub use super::nonlinear_time_of_impact::{
        nonlinear_time_of_impact_ball_ball, nonlinear_time_of_impact_composite_shape_shape,
        nonlinear_time_of_impact_shape_composite_shape,
        nonlinear_time_of_impact_support_map_support_map,
        nonlinear_time_of_impact_support_map_support_map_with_closest_points_function,
    };
    pub use super::point::local_point_projection_on_support_map;
    #[cfg(feature = "dim3")]
    pub use super::ray::local_ray_intersection_with_triangle;
    pub use super::ray::{
        line_toi_with_halfspace, local_ray_intersection_with_support_map_with_params,
        ray_toi_with_ball, ray_toi_with_halfspace,
    };
    pub use super::time_of_impact::{
        time_of_impact_ball_ball, time_of_impact_composite_shape_shape,
        time_of_impact_halfspace_support_map, time_of_impact_shape_composite_shape,
        time_of_impact_support_map_halfspace, time_of_impact_support_map_support_map,
    };
}
