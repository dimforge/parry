use crate::math::{Pose, Real};
use crate::query::{sat, ClosestPoints, PointQuery};
use crate::shape::{Cuboid, SupportMap};

/// Closest points between two cuboids.
#[inline]
pub fn closest_points_cuboid_cuboid(
    pos12: &Pose,
    cuboid1: &Cuboid,
    cuboid2: &Cuboid,
    margin: Real,
) -> ClosestPoints {
    let pos21 = pos12.inverse();

    let sep1 = sat::cuboid_cuboid_find_local_separating_normal_oneway(cuboid1, cuboid2, pos12);
    if sep1.0 > margin {
        return ClosestPoints::Disjoint;
    }

    let sep2 = sat::cuboid_cuboid_find_local_separating_normal_oneway(cuboid2, cuboid1, &pos21);
    if sep2.0 > margin {
        return ClosestPoints::Disjoint;
    }

    #[cfg(feature = "dim2")]
    let sep3 = (-Real::MAX, crate::math::Vector::Y); // This case does not exist in 2D.
    #[cfg(feature = "dim3")]
    let sep3 = sat::cuboid_cuboid_find_local_separating_edge_twoway(cuboid1, cuboid2, pos12);
    if sep3.0 > margin {
        return ClosestPoints::Disjoint;
    }

    if sep1.0 <= 0.0 && sep2.0 <= 0.0 && sep3.0 <= 0.0 {
        return ClosestPoints::Intersecting;
    }

    // The best separating axis is face-vertex.
    if sep1.0 >= sep2.0 && sep1.0 >= sep3.0 {
        // println!("AA: {:?}", sep1);

        let pt2_1 = cuboid2.support_point(pos12, -sep1.1);
        // Use double-projection to converge to the closest points even in ambiguous cases
        // where both cuboids have some parallel features.
        let proj1 = cuboid1.project_local_point(pt2_1, true);
        let proj2 = cuboid2.project_local_point(pos21 * proj1.point, true);
        let pt2_1 = pos12 * proj2.point;
        if (proj1.point - pt2_1).length_squared() > margin * margin {
            return ClosestPoints::Disjoint;
        } else {
            return ClosestPoints::WithinMargin(proj1.point, proj2.point);
        }
    }

    // The best separating axis is vertex-face.
    if sep2.0 >= sep1.0 && sep2.0 >= sep3.0 {
        // println!("BB: {:?}", sep2);

        // Use double-projection to converge to the closest points even in ambiguous cases
        // where both cuboids have some parallel features.
        let pt1_2 = cuboid1.support_point(&pos21, -sep2.1);
        let proj2 = cuboid2.project_local_point(pt1_2, true);
        let proj1 = cuboid1.project_local_point(pos12 * proj2.point, true);
        let pt1_2 = pos21 * proj1.point;

        if (proj2.point - pt1_2).length_squared() > margin * margin {
            return ClosestPoints::Disjoint;
        } else {
            return ClosestPoints::WithinMargin(proj1.point, proj2.point);
        }
    }

    // The best separating axis is edge-edge.
    #[cfg(feature = "dim3")]
    if sep3.0 >= sep2.0 && sep3.0 >= sep1.0 {
        // println!("AA: {:?}, BB: {:?}, CC: {:?}", sep1, sep2, sep3);

        // To compute the actual distance, we need to compute the closest
        // points between the two edges that generated the separating axis.
        let edge1 = cuboid1.local_support_edge_segment(sep3.1);
        let edge2 = cuboid2.local_support_edge_segment(pos21.rotation * -sep3.1);
        let (_, loc2) = super::closest_points_segment_segment_with_locations(pos12, &edge1, &edge2);

        // Use double-projection to converge to the closest points even in ambiguous cases
        // where both cuboids have some parallel features.
        let pt2_1 = pos12 * edge2.point_at(&loc2);
        let proj1 = cuboid1.project_local_point(pt2_1, true);
        let proj2 = cuboid2.project_local_point(pos21 * proj1.point, true);
        let pt2_1 = pos12 * proj2.point;
        if (proj1.point - pt2_1).length_squared() > margin * margin {
            return ClosestPoints::Disjoint;
        } else {
            return ClosestPoints::WithinMargin(proj1.point, proj2.point);
        }
    }

    unreachable!()
}
