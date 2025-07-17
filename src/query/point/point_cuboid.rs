use crate::bounding_volume::Aabb;
use crate::math::{Point, Real};
use crate::query::point::point_query::QueryOptions;
use crate::query::{PointProjection, PointQuery};
use crate::shape::{Cuboid, FeatureId};

impl PointQuery for Cuboid {
    #[inline]
    fn project_local_point(
        &self,
        pt: &Point<Real>,
        solid: bool,
        options: &dyn QueryOptions,
    ) -> PointProjection {
        let dl = Point::from(-self.half_extents);
        let ur = Point::from(self.half_extents);
        Aabb::new(dl, ur).project_local_point(pt, solid, options)
    }

    #[inline]
    fn project_local_point_and_get_feature(
        &self,
        pt: &Point<Real>,
        options: &dyn QueryOptions,
    ) -> (PointProjection, FeatureId) {
        let dl = Point::from(-self.half_extents);
        let ur = Point::from(self.half_extents);
        Aabb::new(dl, ur).project_local_point_and_get_feature(pt, options)
    }

    #[inline]
    fn distance_to_local_point(
        &self,
        pt: &Point<Real>,
        solid: bool,
        options: &dyn QueryOptions,
    ) -> Real {
        let dl = Point::from(-self.half_extents);
        let ur = Point::from(self.half_extents);
        Aabb::new(dl, ur).distance_to_local_point(pt, solid, options)
    }

    #[inline]
    fn contains_local_point(&self, pt: &Point<Real>, _options: &dyn QueryOptions) -> bool {
        let dl = Point::from(-self.half_extents);
        let ur = Point::from(self.half_extents);
        Aabb::new(dl, ur).contains_local_point(pt)
    }
}
