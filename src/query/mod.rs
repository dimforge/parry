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
#[cfg(feature = "std")]
pub use self::contact_manifolds::{
    ContactManifold, ContactManifoldsWorkspace, TrackedContact, TypedWorkspaceData, WorkspaceData,
};
pub use self::default_query_dispatcher::DefaultQueryDispatcher;
pub use self::distance::distance;
pub use self::error::Unsupported;
pub use self::intersection_test::intersection_test;
pub use self::nonlinear_time_of_impact::{nonlinear_time_of_impact, NonlinearRigidMotion};
pub use self::point::{PointProjection, PointQuery, PointQueryWithLocation};
#[cfg(feature = "std")]
pub use self::query_dispatcher::PersistentQueryDispatcher;
pub use self::query_dispatcher::{QueryDispatcher, QueryDispatcherChain};
pub use self::ray::{Ray, RayCast, RayIntersection, SimdRay};
pub use self::time_of_impact::{time_of_impact, TOIStatus, TOI};

mod clip;
pub mod closest_points;
pub mod contact;
#[cfg(feature = "std")]
mod contact_manifolds;
mod default_query_dispatcher;
mod distance;
#[cfg(feature = "std")]
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
#[cfg(feature = "std")]
pub mod visitors;

/// Queries dedicated to specific pairs of shapes.
pub mod details {
    pub use super::clip::*;
    pub use super::closest_points::*;
    pub use super::contact::*;
    #[cfg(feature = "std")]
    pub use super::contact_manifolds::*;
    pub use super::distance::*;
    pub use super::intersection_test::*;
    pub use super::nonlinear_time_of_impact::*;
    pub use super::point::*;
    pub use super::ray::*;
    pub use super::time_of_impact::*;
}
