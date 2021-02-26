use crate::math::{Point, Real, Vector};
use crate::query::{PointProjection, PointQuery};
use crate::shape::{Cone, FeatureId, Segment};
use na;

impl PointQuery for Cone {
    #[inline]
    fn project_local_point(&self, pt: &Point<Real>, solid: bool) -> PointProjection {
        // Project on the basis.
        let mut dir_from_basis_center = pt.coords.xz();
        let planar_dist_from_basis_center = dir_from_basis_center.normalize_mut();

        if planar_dist_from_basis_center <= crate::math::DEFAULT_EPSILON {
            dir_from_basis_center = na::Vector2::x();
        }

        let projection_on_basis = Point::new(pt.coords.x, -self.half_height, pt.coords.z);

        if pt.y < -self.half_height && planar_dist_from_basis_center <= self.radius {
            // The projection is on the basis.
            return PointProjection::new(false, projection_on_basis);
        }

        // Project on the basis circle.
        let proj2d = dir_from_basis_center * self.radius;
        let projection_on_basis_circle = Point::new(proj2d[0], -self.half_height, proj2d[1]);

        // Project on the conic side.
        // TODO: we could solve this in 2D using the plane passing through the cone axis and the conic_side_segment to save some computation.
        let apex_point = Point::new(0.0, self.half_height, 0.0);
        let conic_side_segment = Segment::new(apex_point, projection_on_basis_circle);
        let conic_side_segment_dir = conic_side_segment.scaled_direction();
        let mut proj = conic_side_segment.project_local_point(pt, true);

        let apex_to_basis_center = Vector::new(0.0, -2.0 * self.half_height, 0.0);

        // Now determine if the point is inside of the cone.
        if pt.y >= -self.half_height
            && pt.y <= self.half_height
            && conic_side_segment_dir
                .cross(&(pt - apex_point))
                .dot(&conic_side_segment_dir.cross(&apex_to_basis_center))
                >= 0.0
        {
            if solid {
                PointProjection::new(true, *pt)
            } else {
                // We are inside of the cone, so the correct projection is
                // either on the basis of the cone, or on the conic side.
                if (proj.point - pt).norm_squared() > (projection_on_basis - pt).norm_squared() {
                    PointProjection::new(true, projection_on_basis)
                } else {
                    proj.is_inside = true;
                    proj
                }
            }
        } else {
            // We are outside of the cone, return the computed proj
            // as-is.
            proj
        }
    }

    #[inline]
    fn project_local_point_and_get_feature(
        &self,
        pt: &Point<Real>,
    ) -> (PointProjection, FeatureId) {
        // TODO: get the actual feature.
        (self.project_local_point(pt, false), FeatureId::Unknown)
    }
}
