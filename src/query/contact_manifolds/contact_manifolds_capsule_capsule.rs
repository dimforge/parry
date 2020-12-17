use crate::math::{Isometry, Real, Vector};
use crate::query::{ContactManifold, KinematicsCategory, TrackedContact};
#[cfg(feature = "dim2")]
use crate::shape::CuboidFeature;
#[cfg(feature = "dim2")]
use crate::shape::SegmentPointLocation;
use crate::shape::{Capsule, Shape};
use approx::AbsDiffEq;
use na::Unit;

pub fn contact_manifold_capsule_capsule_shapes<ManifoldData, ContactData>(
    pos12: &Isometry<Real>,
    shape1: &dyn Shape,
    shape2: &dyn Shape,
    prediction: Real,
    manifold: &mut ContactManifold<ManifoldData, ContactData>,
) where
    ContactData: Default + Copy,
{
    if let (Some(capsule1), Some(capsule2)) = (shape1.as_capsule(), shape2.as_capsule()) {
        contact_manifold_capsule_capsule(pos12, capsule1, capsule2, prediction, manifold);
    }
}

#[cfg(feature = "dim2")]
pub fn contact_manifold_capsule_capsule<'a, ManifoldData, ContactData>(
    pos12: &Isometry<Real>,
    capsule1: &'a Capsule,
    capsule2: &'a Capsule,
    prediction: Real,
    manifold: &mut ContactManifold<ManifoldData, ContactData>,
) where
    ContactData: Default + Copy,
{
    // FIXME: the contact kinematics is not correctly set here.
    // We use the common "Point-Plane" kinematics with zero radius everytime.
    // Instead we should select point/point ore point-plane (with non-zero
    // radius for the point) depending on the features involved in the contact.
    let seg1 = capsule1.segment;
    let seg2_1 = capsule2.segment.transformed(&pos12);
    let (loc1, loc2) = crate::query::details::closest_points_segment_segment_with_locations_nD(
        (&seg1.a, &seg1.b),
        (&seg2_1.a, &seg2_1.b),
    );

    // We do this clone to perform contact tracking and transfer impulses.
    // FIXME: find a more efficient way of doing this.
    let old_manifold_points = manifold.points.clone();
    manifold.clear();

    let fid1 = if let SegmentPointLocation::OnVertex(v1) = loc1 {
        v1 as u8 * 2
    } else {
        1
    };
    let fid2 = if let SegmentPointLocation::OnVertex(v2) = loc2 {
        v2 as u8 * 2
    } else {
        1
    };

    let bcoords1 = loc1.barycentric_coordinates();
    let bcoords2 = loc2.barycentric_coordinates();
    let local_p1 = seg1.a * bcoords1[0] + seg1.b.coords * bcoords1[1];
    let local_p2_1 = seg2_1.a * bcoords2[0] + seg2_1.b.coords * bcoords2[1];

    let local_n1 =
        Unit::try_new(local_p2_1 - local_p1, f32::default_epsilon()).unwrap_or(Vector::y_axis());
    let dist = (local_p2_1 - local_p1).dot(&local_n1) - capsule1.radius - capsule2.radius;

    if dist <= prediction {
        let local_n2 = pos12.inverse_transform_unit_vector(&-local_n1);
        let local_p2 = pos12.inverse_transform_point(&local_p2_1);
        let contact = TrackedContact::new(local_p1, local_p2, fid1, fid2, dist);
        manifold.points.push(contact);

        manifold.local_n1 = *local_n1;
        manifold.local_n2 = *local_n2;
        manifold.kinematics.category = KinematicsCategory::PlanePoint;
        manifold.kinematics.radius1 = 0.0;
        manifold.kinematics.radius2 = 0.0;
    } else {
        // No contact within tolerance.
        return;
    }

    if let (Some(dir1), Some(dir2)) = (seg1.direction(), seg2_1.direction()) {
        if dir1.dot(&dir2).abs() >= crate::utils::COS_FRAC_PI_8
            && dir1.dot(&local_n1).abs() < crate::utils::SIN_FRAC_PI_8
        {
            // Capsules axes are almost parallel and are almost perpendicular to the normal.
            // Find a second contact point.
            if let Some((clip_a, clip_b)) = crate::query::details::clip_segment_segment_with_normal(
                (seg1.a, seg1.b),
                (seg2_1.a, seg2_1.b),
                *local_n1,
            ) {
                let contact =
                    if (clip_a.0 - local_p1).norm_squared() > f32::default_epsilon() * 100.0 {
                        // Use clip_a as the second contact.
                        TrackedContact::new(
                            clip_a.0,
                            pos12.inverse_transform_point(&clip_a.1),
                            clip_a.2 as u8,
                            clip_a.3 as u8,
                            (clip_a.1 - clip_a.0).dot(&local_n1),
                        )
                    } else {
                        // Use clip_b as the second contact.
                        TrackedContact::new(
                            clip_b.0,
                            pos12.inverse_transform_point(&clip_b.1),
                            clip_b.2 as u8,
                            clip_b.3 as u8,
                            (clip_b.1 - clip_b.0).dot(&local_n1),
                        )
                    };

                manifold.points.push(contact);
            }
        }
    }

    for point in &mut manifold.points {
        point.local_p1 += manifold.local_n1 * capsule1.radius;
        point.local_p2 += manifold.local_n2 * capsule2.radius;
        point.dist -= capsule1.radius + capsule2.radius;
    }

    manifold.match_contacts(&old_manifold_points);
    manifold.sort_contacts(prediction);
}

#[cfg(feature = "dim3")]
pub fn contact_manifold_capsule_capsule<'a, ManifoldData, ContactData>(
    pos12: &Isometry<Real>,
    capsule1: &'a Capsule,
    capsule2: &'a Capsule,
    prediction: Real,
    manifold: &mut ContactManifold<ManifoldData, ContactData>,
) where
    ContactData: Default + Copy,
{
    let seg1 = capsule1.segment;
    let seg2_1 = capsule2.segment.transformed(&pos12);
    let (loc1, loc2) =
        crate::query::closest_points::closest_points_segment_segment_with_locations_nD(
            (&seg1.a, &seg1.b),
            (&seg2_1.a, &seg2_1.b),
        );

    let bcoords1 = loc1.barycentric_coordinates();
    let bcoords2 = loc2.barycentric_coordinates();
    let local_p1 = seg1.a * bcoords1[0] + seg1.b.coords * bcoords1[1];
    let local_p2_1 = seg2_1.a * bcoords2[0] + seg2_1.b.coords * bcoords2[1];

    let local_n1 =
        Unit::try_new(local_p2_1 - local_p1, f32::default_epsilon()).unwrap_or(Vector::y_axis());
    let dist = (local_p2_1 - local_p1).dot(&local_n1) - capsule1.radius - capsule2.radius;

    if dist <= prediction {
        let local_n2 = pos12.inverse_transform_unit_vector(&-local_n1);
        let contact = TrackedContact::new(
            local_p1 + *local_n1 * capsule1.radius,
            pos12.inverse_transform_point(&local_p2_1) + *local_n2 * capsule2.radius,
            0,
            0,
            dist,
        );

        if manifold.points.len() != 0 {
            manifold.points[0].copy_geometry_from(contact);
        } else {
            manifold.points.push(contact);
        }

        manifold.local_n1 = *local_n1;
        manifold.local_n2 = *local_n2;
        manifold.kinematics.category = KinematicsCategory::PlanePoint;
        manifold.kinematics.radius1 = 0.0;
        manifold.kinematics.radius2 = 0.0;
    } else {
        manifold.clear();
    }

    manifold.sort_contacts(prediction);
}
