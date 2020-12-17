use crate::math::{Point, Real};
use crate::query::{PointProjection, PointQuery};
use crate::shape::{FeatureId, HalfSpace};
use na;

impl PointQuery for HalfSpace {
    #[inline]
    fn project_local_point(&self, pt: &Point<Real>, solid: bool) -> PointProjection {
        let d = self.normal.dot(&pt.coords);
        let inside = d <= na::zero::<Real>();

        if inside && solid {
            PointProjection::new(true, *pt)
        } else {
            PointProjection::new(inside, *pt + (-*self.normal * d))
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
        let dist = self.normal.dot(&pt.coords);

        if dist < na::zero::<Real>() && solid {
            na::zero::<Real>()
        } else {
            // This will automatically be negative if the point is inside.
            dist
        }
    }

    #[inline]
    fn contains_local_point(&self, pt: &Point<Real>) -> bool {
        self.normal.dot(&pt.coords) <= na::zero::<Real>()
    }
}
