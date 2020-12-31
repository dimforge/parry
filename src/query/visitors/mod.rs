//! Visitors for performing geometric queries exploiting spatial partitioning data structures.

pub use self::aabb_sets_interferences_collector::AABBSetsInterferencesCollector;
pub use self::bounding_volume_intersections_visitor::BoundingVolumeIntersectionsVisitor;
pub use self::composite_closest_point_visitor::CompositeClosestPointVisitor;
pub use self::composite_point_containment_test::CompositePointContainmentTest;
pub use self::point_intersections_visitor::PointIntersectionsVisitor;
pub use self::ray_intersections_visitor::RayIntersectionsVisitor;

mod aabb_sets_interferences_collector;
mod bounding_volume_intersections_visitor;
mod composite_closest_point_visitor;
mod composite_point_containment_test;
mod point_intersections_visitor;
mod ray_intersections_visitor;
