use crate::math::{Pose, Real, Vector};
use crate::query::contact_manifolds::{NormalConstraints, NormalConstraintsPair};
use crate::query::{
    self,
    gjk::{GJKResult, VoronoiSimplex},
    ContactManifold, TrackedContact,
};
use crate::shape::{PackedFeatureId, PolygonalFeature, PolygonalFeatureMap, Shape};

/// Computes the contact manifold between two convex shapes implementing the `PolygonalSupportMap`
/// trait, both represented as `Shape` trait-objects.
pub fn contact_manifold_pfm_pfm_shapes<ManifoldData, ContactData>(
    pos12: &Pose,
    shape1: &dyn Shape,
    shape2: &dyn Shape,
    normal_constraints1: Option<&dyn NormalConstraints>,
    normal_constraints2: Option<&dyn NormalConstraints>,
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
            normal_constraints1,
            pfm2,
            border_radius2,
            normal_constraints2,
            prediction,
            manifold,
        );
    }
}

/// Computes the contact manifold between two convex shapes implementing the `PolygonalSupportMap` trait.
pub fn contact_manifold_pfm_pfm<'a, ManifoldData, ContactData, S1, S2>(
    pos12: &Pose,
    pfm1: &'a S1,
    border_radius1: Real,
    normal_constraints1: Option<&dyn NormalConstraints>,
    pfm2: &'a S2,
    border_radius2: Real,
    normal_constraints2: Option<&dyn NormalConstraints>,
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
    // Never for a curved surface: when it rolls, its material contact point barely moves
    // relative to the other shape (so the update would pass) but leaves the true contact point,
    // up to `r·sin(1°)` behind it, and that lever arm makes the normal force spin the roll up.
    let curved = is_curved(pfm1, border_radius1) || is_curved(pfm2, border_radius2);
    if !curved && manifold.try_update_contacts_eps(pos12, crate::utils::COS_1_DEGREES, 1.0e-6) {
        return;
    }

    let init_dir = (manifold.local_n1).try_normalize();
    let total_prediction = prediction + border_radius1 + border_radius2;
    let contact = query::details::contact_support_map_support_map_with_params(
        pos12,
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
            #[cfg(feature = "dim3")]
            let (p1, p2_1, dir) = if curved {
                refine_normal(pos12, pfm1, pfm2, init_dir, (p1, p2_1, dir))
            } else {
                (p1, p2_1, dir)
            };

            let mut local_n1 = dir;
            let mut local_n2 = pos12.rotation.inverse() * -dir;
            let dist = (p2_1 - p1).dot(local_n1);

            if !(normal_constraints1, normal_constraints2).project_local_normals(
                pos12,
                &mut local_n1,
                &mut local_n2,
            ) {
                // The contact got completely discarded by the normal correction.
                return;
            }

            let mut feature1 = PolygonalFeature::default();
            let mut feature2 = PolygonalFeature::default();
            pfm1.local_support_feature_toward(local_n1, p1, &mut feature1);
            pfm2.local_support_feature_toward(
                local_n2,
                pos12.inverse_transform_point(p2_1),
                &mut feature2,
            );

            PolygonalFeature::contacts(
                pos12,
                &pos12.inverse(),
                local_n1,
                local_n2,
                &feature1,
                &feature2,
                manifold,
                false,
            );

            if (cfg!(feature = "dim3") || cfg!(feature = "dim2") && manifold.points.is_empty())
                // If normal constraints changed the GJK direction, it is no longer valid so we cant use it for this additional contact.
                && local_n1 == dir
            {
                let mut contact = TrackedContact::new(
                    p1,
                    pos12.inverse_transform_point(p2_1),
                    PackedFeatureId::UNKNOWN, // TODO: We don't know what features are involved.
                    PackedFeatureId::UNKNOWN,
                    (p2_1 - p1).dot(local_n1),
                );

                // On a curved surface the GJK witness point is only near the contact (a few
                // tenths of a millimetre to millimetres off it, across the contact line of a
                // cylinder's side), which gives its normal force a lever arm about the axis.
                // When the features gave a contact line, the extra point is moved onto it,
                // unless that line is farther than the witness (a capsule's segment on the far
                // side of a box edge), which would lose the actual contact.
                if curved && manifold.points.len() == 2 {
                    let (a, b) = (&manifold.points[0], &manifold.points[1]);
                    let ab = b.local_p1 - a.local_p1;
                    let t = if ab.length_squared() > 0.0 {
                        ((p1 - a.local_p1).dot(ab) / ab.length_squared()).clamp(0.0, 1.0)
                    } else {
                        0.5
                    };
                    let dist = a.dist + (b.dist - a.dist) * t;
                    if dist <= contact.dist + 1.0e-4 * (1.0 + ab.length()) {
                        contact.local_p1 = a.local_p1 + ab * t;
                        contact.local_p2 = a.local_p2 + (b.local_p2 - a.local_p2) * t;
                        contact.dist = dist;
                    }
                }

                manifold.points.push(contact);
            }

            if normal_constraints1.is_some() || normal_constraints2.is_some() {
                // HACK: some normal correction can lead to very incorrect penetration
                //       depth, e.g., if the other object extends very far toward that direction.
                //       This is caused by the locality of the convex/convex check.
                //       I haven’t found a good mathematically robust approach to account for
                //       that locally, so for now, we eliminate points that are large divergence
                //       relative to the unconstrained penetration distance.
                manifold
                    .points
                    .retain(|pt| dist >= 0.0 || pt.dist >= 0.0 || pt.dist >= dist * 5.0);
            }

            // Adjust points to take the radius into account.
            if border_radius1 != 0.0 || border_radius2 != 0.0 {
                for contact in &mut manifold.points {
                    contact.local_p1 += local_n1 * border_radius1;
                    contact.local_p2 += local_n2 * border_radius2;
                    contact.dist -= border_radius1 + border_radius2;
                }
            }

            manifold.local_n1 = local_n1;
            manifold.local_n2 = local_n2;
        }
        GJKResult::NoIntersection(dir) => {
            // Use the manifold normal as a cache.
            manifold.local_n1 = dir;
        }
        _ => {
            // Reset the cached direction.
            manifold.local_n1 = Vector::ZERO;
        }
    }

    // Transfer impulses.
    manifold.match_contacts(&old_manifold_points);
}

/// A curved surface: a cylinder's or cone's side, or any rounded shape.
fn is_curved<S: ?Sized + PolygonalFeatureMap>(pfm: &S, border_radius: Real) -> bool {
    border_radius != 0.0 || pfm.has_curved_surface()
}

/// GJK's and EPA's normal for a curved shape can be tens of degrees off when the shapes are
/// barely apart or barely overlap. The face either shape shows toward that normal is a
/// separating-axis candidate, and if the shapes are farther apart along it, it's the better
/// normal (always so for a curved shape resting on a flat face).
#[cfg(feature = "dim3")]
fn refine_normal<S1, S2>(
    pos12: &Pose,
    pfm1: &S1,
    pfm2: &S2,
    prev_dir: Option<Vector>,
    (p1, p2_1, dir): (Vector, Vector, Vector),
) -> (Vector, Vector, Vector)
where
    S1: ?Sized + PolygonalFeatureMap,
    S2: ?Sized + PolygonalFeatureMap,
{
    let face_normal = |f: &PolygonalFeature| {
        (f.num_vertices >= 3)
            .then(|| (f.vertices[1] - f.vertices[0]).cross(f.vertices[2] - f.vertices[0]))
            .and_then(|n| n.try_normalize())
    };
    let faces_toward = |n: Vector| {
        let mut feature1 = PolygonalFeature::default();
        let mut feature2 = PolygonalFeature::default();
        pfm1.local_support_feature(n, &mut feature1);
        pfm2.local_support_feature(pos12.rotation.inverse() * -n, &mut feature2);
        [
            face_normal(&feature1),
            face_normal(&feature2).map(|n2| pos12.rotation * -n2),
        ]
    };
    // The faces toward GJK's normal, toward the other shape's center, and last step's normal
    // with the faces toward it.
    let [a, b] = faces_toward(dir);
    let [c, d] = pos12
        .translation
        .try_normalize()
        .map_or([None, None], faces_toward);
    let [e, f] = prev_dir.map_or([None, None], faces_toward);
    let candidates = [a, b, c, d, prev_dir, e, f];

    // Separation along `n` measured on the shapes themselves, not GJK's witness points.
    let separation = |n: Vector| {
        let a = pfm1.local_support_point(n);
        let b = pos12 * pfm2.local_support_point(pos12.rotation.inverse() * -n);
        ((b - a).dot(n), b)
    };
    let mut best = (separation(dir).0, (p1, p2_1, dir));
    // Each axis is tried both ways since GJK's normal may be too far off to orient it.
    for n in candidates.into_iter().flatten().flat_map(|n| [n, -n]) {
        let (sep, b) = separation(n);
        if sep > best.0 + 1.0e-6 {
            best = (sep, (b - n * sep, b, n));
        }
    }
    best.1
}
