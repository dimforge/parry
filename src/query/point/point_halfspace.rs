use crate::math::*;
use crate::query::{PointProjection, PointQuery};
use crate::shape::{FeatureId, HalfSpace};

impl PointQuery for HalfSpace {
    #[inline]
    fn project_local_point(&self, pt: &Point, solid: bool) -> PointProjection {
        let d = self.normal.dot(pt.as_vector());
        let inside = d <= 0.0;

        if inside && solid {
            PointProjection::new(true, *pt)
        } else {
            PointProjection::new(inside, *pt + (-self.normal.into_inner() * d))
        }
    }

    #[inline]
    fn project_local_point_and_get_feature(&self, pt: &Point) -> (PointProjection, FeatureId) {
        (self.project_local_point(pt, false), FeatureId::Face(0))
    }

    #[inline]
    fn distance_to_local_point(&self, pt: &Point, solid: bool) -> Real {
        let dist = self.normal.dot(pt.as_vector());

        if dist < 0.0 && solid {
            0.0
        } else {
            // This will automatically be negative if the point is inside.
            dist
        }
    }

    #[inline]
    fn contains_local_point(&self, pt: &Point) -> bool {
        self.normal.dot(pt.as_vector()) <= 0.0
    }
}
