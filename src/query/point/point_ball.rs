use na::{self, ComplexField};

use crate::math::{Point, Real};
use crate::query::{PointProjection, PointQuery};
use crate::shape::{Ball, FeatureId};

impl PointQuery for Ball {
    #[inline]
    fn project_local_point(&self, pt: &Point<Real>, solid: bool) -> PointProjection {
        let distance_squared = pt.coords.norm_squared();

        let inside = distance_squared <= self.radius * self.radius;

        if inside && solid {
            PointProjection::new(true, *pt)
        } else {
            let proj =
                Point::from(pt.coords * (self.radius / ComplexField::sqrt(distance_squared)));
            PointProjection::new(inside, proj)
        }
    }

    #[inline]
    fn project_local_point_and_get_feature(
        &self,
        pt: &Point<Real>,
    ) -> (PointProjection, FeatureId) {
        (self.project_local_point(pt, false), FeatureId::Face(0))
    }

    #[inline]
    fn distance_to_local_point(&self, pt: &Point<Real>, solid: bool) -> Real {
        let dist = pt.coords.norm() - self.radius;

        if solid && dist < 0.0 {
            0.0
        } else {
            dist
        }
    }

    #[inline]
    fn contains_local_point(&self, pt: &Point<Real>) -> bool {
        pt.coords.norm_squared() <= self.radius * self.radius
    }
}
