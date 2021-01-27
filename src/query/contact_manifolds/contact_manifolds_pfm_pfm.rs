use crate::math::{Isometry, Real};
use crate::query::{
    self,
    gjk::{GJKResult, VoronoiSimplex},
    ContactManifold, TrackedContact,
};
use crate::shape::{PolygonalFeature, PolygonalFeatureMap, Shape};
use na::Unit;

/// Computes the contact manifold between two convex shapes implementing the `PolygonalSupportMap`
/// trait, both represented as `Shape` trait-objects.
pub fn contact_manifold_pfm_pfm_shapes<ManifoldData, ContactData>(
    pos12: &Isometry<Real>,
    shape1: &dyn Shape,
    shape2: &dyn Shape,
    prediction: Real,
    manifold: &mut ContactManifold<ManifoldData, ContactData>,
) where
    ManifoldData: Default,
    ContactData: Default + Copy,
{
    if let (Some((pfm1, border_radius1)), Some((pfm2, border_radius2))) = (
        shape1.as_polygonal_feature_map(),
        shape2.as_polygonal_feature_map(),
    ) {
        contact_manifold_pfm_pfm(
            pos12,
            pfm1,
            border_radius1,
            pfm2,
            border_radius2,
            prediction,
            manifold,
        );
    }
}

/// Computes the contact manifold between two convex shapes implementing the `PolygonalSupportMap` trait.
pub fn contact_manifold_pfm_pfm<'a, ManifoldData, ContactData, S1, S2>(
    pos12: &Isometry<Real>,
    pfm1: &'a S1,
    border_radius1: Real,
    pfm2: &'a S2,
    border_radius2: Real,
    prediction: Real,
    manifold: &mut ContactManifold<ManifoldData, ContactData>,
) where
    S1: ?Sized + PolygonalFeatureMap,
    S2: ?Sized + PolygonalFeatureMap,
    ManifoldData: Default,
    ContactData: Default + Copy,
{
    // We use very small thresholds for the manifold update because something to high would
    // cause numerical drifts with the effect of introducing bumps in
    // what should have been smooth rolling motions.
    if manifold.try_update_contacts_eps(&pos12, crate::utils::COS_1_DEGREES, 1.0e-6) {
        return;
    }

    let init_dir = Unit::try_new(manifold.local_n1, crate::math::DEFAULT_EPSILON);
    let total_prediction = prediction + border_radius1 + border_radius2;
    let contact = query::details::contact_support_map_support_map_with_params(
        &pos12,
        pfm1,
        pfm2,
        total_prediction,
        &mut VoronoiSimplex::new(),
        init_dir,
    );

    let old_manifold_points = manifold.points.clone();
    manifold.clear();

    match contact {
        GJKResult::ClosestPoints(p1, p2_1, dir) => {
            let local_n1 = dir;
            let local_n2 = pos12.inverse_transform_unit_vector(&-dir);
            let mut feature1 = PolygonalFeature::default();
            let mut feature2 = PolygonalFeature::default();
            pfm1.local_support_feature(&local_n1, &mut feature1);
            pfm2.local_support_feature(&local_n2, &mut feature2);

            PolygonalFeature::contacts(
                pos12,
                &pos12.inverse(),
                &local_n1,
                &local_n2,
                &feature1,
                &feature2,
                total_prediction,
                manifold,
                false,
            );

            if cfg!(feature = "dim3") || (cfg!(feature = "dim2") && manifold.points.is_empty()) {
                let contact = TrackedContact::new(
                    p1,
                    pos12.inverse_transform_point(&p2_1),
                    u32::MAX, // We don't know what features are involved.
                    u32::MAX,
                    (p2_1 - p1).dot(&dir),
                );
                manifold.points.push(contact);
            }

            // Adjust points to take the radius into account.
            if border_radius1 != 0.0 || border_radius2 != 0.0 {
                for contact in &mut manifold.points {
                    contact.local_p1 += *local_n1 * border_radius1;
                    contact.local_p2 += *local_n2 * border_radius2;
                    contact.dist -= border_radius1 + border_radius2;
                }
            }

            manifold.local_n1 = *local_n1;
            manifold.local_n2 = *local_n2;
        }
        GJKResult::NoIntersection(dir) => {
            // Use the manifold normal as a cache.
            manifold.local_n1 = *dir;
        }
        _ => {
            // Reset the cached direction.
            manifold.local_n1.fill(0.0);
        }
    }

    // Transfer impulses.
    manifold.match_contacts(&old_manifold_points);
}
