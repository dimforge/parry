use crate::math::*;
use crate::query::{PointProjection, PointQuery};
use crate::shape::{Cylinder, FeatureId};

#[cfg(feature = "linalg-glam")]
use glam::Vec3Swizzles;

impl PointQuery for Cylinder {
    #[inline]
    fn project_local_point(&self, pt: &Point, solid: bool) -> PointProjection {
        // Project on the basis.
        let mut dir_from_basis_center = pt.as_vector().xz();
        let planar_dist_from_basis_center = dir_from_basis_center.normalize_mut();

        if planar_dist_from_basis_center <= DEFAULT_EPSILON {
            dir_from_basis_center = Vector2::x();
        }

        let proj2d = dir_from_basis_center * self.radius;

        if pt.y >= -self.half_height
            && pt.y <= self.half_height
            && planar_dist_from_basis_center <= self.radius
        {
            // The point is inside of the cylinder.
            if solid {
                PointProjection::new(true, *pt)
            } else {
                let dist_to_top = self.half_height - pt.as_vector().y;
                let dist_to_bottom = pt.as_vector().y - (-self.half_height);
                let dist_to_side = self.radius - planar_dist_from_basis_center;

                if dist_to_top < dist_to_bottom && dist_to_top < dist_to_side {
                    let projection_on_top =
                        Point::new(pt.as_vector().x, self.half_height, pt.as_vector().z);
                    PointProjection::new(true, projection_on_top)
                } else if dist_to_bottom < dist_to_top && dist_to_bottom < dist_to_side {
                    let projection_on_bottom =
                        Point::new(pt.as_vector().x, -self.half_height, pt.as_vector().z);
                    PointProjection::new(true, projection_on_bottom)
                } else {
                    let projection_on_side = Point::new(proj2d[0], pt.y, proj2d[1]);
                    PointProjection::new(true, projection_on_side)
                }
            }
        } else {
            // The point is outside of the cylinder.
            if pt.y > self.half_height {
                if planar_dist_from_basis_center <= self.radius {
                    let projection_on_top =
                        Point::new(pt.as_vector().x, self.half_height, pt.as_vector().z);
                    PointProjection::new(false, projection_on_top)
                } else {
                    let projection_on_top_circle =
                        Point::new(proj2d[0], self.half_height, proj2d[1]);
                    PointProjection::new(false, projection_on_top_circle)
                }
            } else if pt.y < -self.half_height {
                // Project on the bottom plane or the bottom circle.
                if planar_dist_from_basis_center <= self.radius {
                    let projection_on_bottom =
                        Point::new(pt.as_vector().x, -self.half_height, pt.as_vector().z);
                    return PointProjection::new(false, projection_on_bottom);
                } else {
                    let projection_on_bottom_circle =
                        Point::new(proj2d[0], -self.half_height, proj2d[1]);
                    return PointProjection::new(false, projection_on_bottom_circle);
                }
            } else {
                // Project on the side.
                let projection_on_side = Point::new(proj2d[0], pt.y, proj2d[1]);
                return PointProjection::new(false, projection_on_side);
            }
        }
    }

    #[inline]
    fn project_local_point_and_get_feature(&self, pt: &Point) -> (PointProjection, FeatureId) {
        // TODO: get the actual feature.
        (self.project_local_point(pt, false), FeatureId::Unknown)
    }
}
