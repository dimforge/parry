use na::{self, ComplexField};

use crate::math::*;
use crate::query::{PointProjection, PointQuery};
use crate::shape::{Ball, FeatureId};

impl PointQuery for Ball {
    #[inline]
    fn project_local_point(&self, pt: &Point, solid: bool) -> PointProjection {
        let distance_squared = pt.as_vector().norm_squared();

        let inside = distance_squared <= self.radius * self.radius;

        if inside && solid {
            PointProjection::new(true, *pt)
        } else {
            let proj =
                Point::from(pt.as_vector() * (self.radius / ComplexField::sqrt(distance_squared)));
            PointProjection::new(inside, proj)
        }
    }

    #[inline]
    fn project_local_point_and_get_feature(&self, pt: &Point) -> (PointProjection, FeatureId) {
        (self.project_local_point(pt, false), FeatureId::Face(0))
    }

    #[inline]
    fn distance_to_local_point(&self, pt: &Point, solid: bool) -> Real {
        let dist = pt.as_vector().norm() - self.radius;

        if solid && dist < 0.0 {
            0.0
        } else {
            dist
        }
    }

    #[inline]
    fn contains_local_point(&self, pt: &Point) -> bool {
        pt.as_vector().norm_squared() <= self.radius * self.radius
    }
}
