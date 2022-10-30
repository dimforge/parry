//! Visitors for performing geometric queries exploiting spatial partitioning data structures.

#[cfg(feature = "std")]
pub use self::aabb_sets_interferences_collector::AabbSetsInterferencesCollector;
#[cfg(feature = "std")]
pub use self::bounding_volume_intersections_simultaneous_visitor::BoundingVolumeIntersectionsSimultaneousVisitor;
#[cfg(feature = "std")]
pub use self::bounding_volume_intersections_visitor::BoundingVolumeIntersectionsVisitor;
#[cfg(feature = "std")]
pub use self::composite_closest_point_visitor::CompositeClosestPointVisitor;
pub use self::composite_point_containment_test::CompositePointContainmentTest;
#[cfg(feature = "std")]
pub use self::point_intersections_visitor::PointIntersectionsVisitor;
#[cfg(feature = "std")]
pub use self::ray_intersections_visitor::RayIntersectionsVisitor;

#[cfg(feature = "std")]
mod aabb_sets_interferences_collector;
#[cfg(feature = "std")]
mod bounding_volume_intersections_simultaneous_visitor;
#[cfg(feature = "std")]
mod bounding_volume_intersections_visitor;
#[cfg(feature = "std")]
mod composite_closest_point_visitor;
mod composite_point_containment_test;
#[cfg(feature = "std")]
mod point_intersections_visitor;
#[cfg(feature = "std")]
mod ray_intersections_visitor;
