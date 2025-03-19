use crate::math::{Isometry, Real};
use crate::query::{sat, Contact, PointQuery};
use crate::shape::{Cuboid, SupportMap};
use approx::AbsDiffEq;
use na::Unit;

/// Contact between two cuboids.
#[inline]
pub fn contact_cuboid_cuboid(
    pos12: &Isometry,
    cuboid1: &Cuboid,
    cuboid2: &Cuboid,
    prediction: Real,
) -> Option<Contact> {
    let pos21 = pos12.inverse();

    let sep1 = sat::cuboid_cuboid_find_local_separating_normal_oneway(cuboid1, cuboid2, pos12);
    if sep1.0 > prediction {
        return None;
    }

    let sep2 = sat::cuboid_cuboid_find_local_separating_normal_oneway(cuboid2, cuboid1, &pos21);
    if sep2.0 > prediction {
        return None;
    }

    #[cfg(feature = "dim2")]
    let sep3 = (-Real::MAX, crate::math::Vector::y()); // This case does not exist in 2D.
    #[cfg(feature = "dim3")]
    let sep3 = sat::cuboid_cuboid_find_local_separating_edge_twoway(cuboid1, cuboid2, pos12);
    if sep3.0 > prediction {
        return None;
    }

    // The best separating axis is face-vertex.
    if sep1.0 >= sep2.0 && sep1.0 >= sep3.0 {
        // To compute the closest points, we need to project the support point
        // from cuboid2 on the support-face of cuboid1. For simplicity, we just
        // project the support point from cuboid2 on cuboid1 itself (not just the face).
        let pt2_1 = cuboid2.support_point(pos12, &-sep1.1);
        let proj1 = cuboid1.project_local_point(&pt2_1, false);

        let separation = (pt2_1 - proj1.point).dot(&sep1.1);
        let normalized_dir = Unit::try_new_and_get(pt2_1 - proj1.point, Real::default_epsilon());
        let normal1;
        let dist;

        // NOTE: we had to recompute the normal because we can't use
        // the separation vector for the case where we have a vertex-vertex contact.
        if separation < 0.0 || normalized_dir.is_none() {
            // Penetration or contact lying on the boundary exactly.
            normal1 = Unit::new_unchecked(sep1.1);
            dist = separation;
        } else {
            let (dir, norm) = normalized_dir.unwrap();
            // No penetration.
            normal1 = dir;
            dist = norm;
        }

        if dist > prediction {
            return None;
        }

        return Some(Contact::new(
            proj1.point,
            pos12.inverse_transform_point(&pt2_1),
            normal1,
            pos12.inverse_transform_unit_vector(&-normal1),
            dist,
        ));
    }

    // The best separating axis is vertex-face.
    if sep2.0 >= sep1.0 && sep2.0 >= sep3.0 {
        // To compute the actual closest points, we need to project the support point
        // from cuboid1 on the support-face of cuboid2. For simplicity, we just
        // project the support point from cuboid1 on cuboid2 itself (not just the face).
        let pt1_2 = cuboid1.support_point(&pos21, &-sep2.1);
        let proj2 = cuboid2.project_local_point(&pt1_2, false);

        let separation = (pt1_2 - proj2.point).dot(&sep2.1);
        let normalized_dir = Unit::try_new_and_get(pt1_2 - proj2.point, Real::default_epsilon());
        let normal2;
        let dist;

        // NOTE: we had to recompute the normal because we can't use
        // the separation vector for the case where we have a vertex-vertex contact.
        if separation < 0.0 || normalized_dir.is_none() {
            // Penetration or contact lying on the boundary exactly.
            normal2 = Unit::new_unchecked(sep2.1);
            dist = separation;
        } else {
            // No penetration.
            let (dir, norm) = normalized_dir.unwrap();
            normal2 = dir;
            dist = norm;
        }

        if dist > prediction {
            return None;
        }

        return Some(Contact::new(
            pos12.transform_point(&pt1_2),
            proj2.point,
            pos12 * -normal2,
            normal2,
            dist,
        ));
    }

    // The best separating axis is edge-edge.
    #[cfg(feature = "dim3")]
    if sep3.0 >= sep2.0 && sep3.0 >= sep1.0 {
        use crate::query::{details, ClosestPoints};
        // To compute the actual distance, we need to compute the closest
        // points between the two edges that generated the separating axis.
        let edge1 = cuboid1.local_support_edge_segment(sep3.1);
        let edge2 = cuboid2.local_support_edge_segment(pos21 * -sep3.1);

        match details::closest_points_segment_segment(pos12, &edge1, &edge2, prediction) {
            ClosestPoints::Disjoint => return None,
            ClosestPoints::WithinMargin(a, b) => {
                let normal1 = Unit::new_unchecked(sep3.1);
                let normal2 = pos12.inverse_transform_unit_vector(&-normal1);
                return Some(Contact::new(a, b, normal1, normal2, sep3.0));
            }
            ClosestPoints::Intersecting => unreachable!(),
        }
    }

    unreachable!()
}
