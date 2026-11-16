use crate::math::{Vector, Vector2};
use crate::query::{PointProjection, PointQuery};
use crate::shape::{Cylinder, FeatureId};

impl PointQuery for Cylinder {
    #[inline]
    fn project_local_point(&self, pt: Vector, solid: bool) -> PointProjection {
        // Project on the basis.
        let (mut dir_from_basis_center, planar_dist_from_basis_center) =
            Vector2::new(pt.x, pt.z).normalize_and_length();

        if planar_dist_from_basis_center <= crate::math::DEFAULT_EPSILON {
            {
                dir_from_basis_center = Vector2::X;
            }
        }

        let proj2d = dir_from_basis_center * self.radius;

        if pt.y >= -self.half_height
            && pt.y <= self.half_height
            && planar_dist_from_basis_center <= self.radius
        {
            // The point is inside of the cylinder.
            if solid {
                PointProjection::new(true, pt)
            } else {
                let dist_to_top = self.half_height - pt.y;
                let dist_to_bottom = pt.y - (-self.half_height);
                let dist_to_side = self.radius - planar_dist_from_basis_center;

                if dist_to_top <= dist_to_bottom && dist_to_top < dist_to_side {
                    let projection_on_top = Vector::new(pt.x, self.half_height, pt.z);
                    PointProjection::new(true, projection_on_top)
                } else if dist_to_bottom < dist_to_top && dist_to_bottom < dist_to_side {
                    let projection_on_bottom = Vector::new(pt.x, -self.half_height, pt.z);
                    PointProjection::new(true, projection_on_bottom)
                } else {
                    let projection_on_side = Vector::new(proj2d[0], pt.y, proj2d[1]);
                    PointProjection::new(true, projection_on_side)
                }
            }
        } else {
            // The point is outside of the cylinder.
            if pt.y > self.half_height {
                if planar_dist_from_basis_center <= self.radius {
                    let projection_on_top = Vector::new(pt.x, self.half_height, pt.z);
                    PointProjection::new(false, projection_on_top)
                } else {
                    let projection_on_top_circle =
                        Vector::new(proj2d[0], self.half_height, proj2d[1]);
                    PointProjection::new(false, projection_on_top_circle)
                }
            } else if pt.y < -self.half_height {
                // Project on the bottom plane or the bottom circle.
                if planar_dist_from_basis_center <= self.radius {
                    let projection_on_bottom = Vector::new(pt.x, -self.half_height, pt.z);
                    PointProjection::new(false, projection_on_bottom)
                } else {
                    let projection_on_bottom_circle =
                        Vector::new(proj2d[0], -self.half_height, proj2d[1]);
                    PointProjection::new(false, projection_on_bottom_circle)
                }
            } else {
                // Project on the side.
                let projection_on_side = Vector::new(proj2d[0], pt.y, proj2d[1]);
                PointProjection::new(false, projection_on_side)
            }
        }
    }

    #[inline]
    fn project_local_point_and_get_feature(&self, pt: Vector) -> (PointProjection, FeatureId) {
        // TODO: get the actual feature.
        (self.project_local_point(pt, false), FeatureId::Unknown)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::Pose;

    #[test]
    fn equal_nearer_caps_project_to_the_boundary() {
        let cylinder = Cylinder::new(0.5, 2.0);
        for point in [Vector::ZERO, Vector::new(0.25, 0.0, 0.25)] {
            let projection = cylinder.project_local_point(point, false);
            assert!(projection.is_inside);
            assert_eq!(projection.point.x, point.x);
            assert_eq!(projection.point.z, point.z);
            assert_eq!(projection.point.y.abs(), 0.5);
            assert_eq!((projection.point - point).length(), 0.5);
            assert_eq!(cylinder.distance_to_local_point(point, false), -0.5);

            let bounded = cylinder
                .project_local_point_with_max_dist(point, false, 0.5)
                .unwrap();
            assert!(bounded.is_inside);
            assert_eq!(bounded.point, projection.point);
            assert!(cylinder
                .project_local_point_with_max_dist(point, false, 0.25)
                .is_none());

            let solid = cylinder.project_local_point(point, true);
            assert!(solid.is_inside);
            assert_eq!(solid.point, point);
            assert_eq!(cylinder.distance_to_local_point(point, true), 0.0);
        }
    }

    #[test]
    fn equal_nearer_caps_project_in_world_space() {
        let cylinder = Cylinder::new(0.5, 2.0);
        let pose = Pose::new(Vector::new(3.0, 4.0, 5.0), Vector::new(0.0, 0.0, 0.7));
        let point = pose.transform_point(Vector::ZERO);
        let projection = cylinder.project_point(&pose, point, false);
        assert!(projection.is_inside);
        assert!(((projection.point - point).length() - 0.5).abs() < 1.0e-5);
        assert!((cylinder.distance_to_point(&pose, point, false) + 0.5).abs() < 1.0e-5);
        let bounded = cylinder
            .project_point_with_max_dist(&pose, point, false, 0.75)
            .unwrap();
        assert!(bounded.is_inside);
        assert_eq!(bounded.point, projection.point);
        assert!(cylinder
            .project_point_with_max_dist(&pose, point, false, 0.25)
            .is_none());
    }

    #[cfg(feature = "alloc")]
    #[test]
    fn equal_nearer_caps_project_through_shared_shape() {
        let cylinder = crate::shape::SharedShape::cylinder(0.5, 2.0);
        let projection = cylinder.project_local_point(Vector::ZERO, false);
        assert!(projection.is_inside);
        assert_eq!(projection.point.length(), 0.5);
        assert_eq!(cylinder.distance_to_local_point(Vector::ZERO, false), -0.5);
        assert!(cylinder
            .project_local_point_with_max_dist(Vector::ZERO, false, 0.75)
            .is_some());
    }

    #[test]
    fn other_nearest_surfaces_are_preserved() {
        let cylinder = Cylinder::new(0.5, 2.0);
        for (point, expected, is_inside) in [
            (Vector::Y * 0.25, Vector::Y * 0.5, true),
            (-Vector::Y * 0.25, -Vector::Y * 0.5, true),
            (Vector::Y, Vector::Y * 0.5, false),
            (Vector::X * 3.0, Vector::X * 2.0, false),
        ] {
            let projection = cylinder.project_local_point(point, false);
            assert_eq!(projection.is_inside, is_inside);
            assert_eq!(projection.point, expected);
        }
    }

    #[test]
    fn nearer_sides_and_side_ties_are_preserved() {
        for (half_height, radius, point) in [
            (2.0, 0.5, Vector::ZERO),
            (0.5, 0.5, Vector::ZERO),
            (0.5, 0.25, Vector::Y * 0.25),
            (0.5, 0.25, -Vector::Y * 0.25),
        ] {
            let projection = Cylinder::new(half_height, radius).project_local_point(point, false);
            assert!(projection.is_inside);
            assert_eq!(projection.point.y, point.y);
            assert_eq!(
                Vector2::new(projection.point.x, projection.point.z).length(),
                radius
            );
            assert_eq!((projection.point - point).length(), radius);
        }
    }
}
