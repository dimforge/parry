#[cfg(feature = "dim2")]
use crate::math::Vector;
use crate::math::{Isometry, Real};
use crate::query::{sat, ContactManifold, KinematicsCategory};
use crate::shape::{Cuboid, CuboidFeature, Shape};

pub fn contact_manifold_cuboid_cuboid_shapes<ManifoldData, ContactData: Default + Copy>(
    pos12: &Isometry<Real>,
    g1: &dyn Shape,
    g2: &dyn Shape,
    prediction: Real,
    manifold: &mut ContactManifold<ManifoldData, ContactData>,
) {
    if let (Some(cuboid1), Some(cuboid2)) = (g1.as_cuboid(), g2.as_cuboid()) {
        contact_manifold_cuboid_cuboid(pos12, cuboid1, cuboid2, prediction, manifold);
    }
}

pub fn contact_manifold_cuboid_cuboid<'a, ManifoldData, ContactData: Default + Copy>(
    pos12: &Isometry<Real>,
    cuboid1: &'a Cuboid,
    cuboid2: &'a Cuboid,
    prediction: Real,
    manifold: &mut ContactManifold<ManifoldData, ContactData>,
) {
    if manifold.try_update_contacts(&pos12) {
        manifold.sort_contacts(prediction);
        return;
    }

    let pos21 = &pos12.inverse();
    let pos12 = &*pos12;

    /*
     *
     * Point-Face
     *
     */
    let sep1 = sat::cuboid_cuboid_find_local_separating_normal_oneway(cuboid1, cuboid2, &pos12);
    if sep1.0 > prediction {
        manifold.clear();
        return;
    }

    let sep2 = sat::cuboid_cuboid_find_local_separating_normal_oneway(cuboid2, cuboid1, &pos21);
    if sep2.0 > prediction {
        manifold.clear();
        return;
    }

    /*
     *
     * Edge-Edge cases
     *
     */
    #[cfg(feature = "dim2")]
    let sep3 = (-Real::MAX, Vector::x()); // This case does not exist in 2D.
    #[cfg(feature = "dim3")]
    let sep3 = sat::cuboid_cuboid_find_local_separating_edge_twoway(cuboid1, cuboid2, &pos12);
    if sep3.0 > prediction {
        manifold.clear();
        return;
    }

    /*
     *
     * Select the best combination of features
     * and get the polygons to clip.
     *
     */
    let mut best_sep = sep1;

    if sep2.0 > sep1.0 && sep2.0 > sep3.0 {
        best_sep = (sep2.0, pos12 * -sep2.1);
    } else if sep3.0 > sep1.0 {
        best_sep = sep3;
    }

    // We do this clone to perform contact tracking and transfer impulses.
    // FIXME: find a more efficient way of doing this.
    let old_manifold_points = manifold.points.clone();
    manifold.clear();
    let local_n2 = pos21 * -best_sep.1;

    // Now the reference feature is from `cuboid1` and the best separation is `best_sep`.
    // Everything must be expressed in the local-space of `cuboid1` for contact clipping.
    let feature1 = cuboid1.support_feature(best_sep.1);
    let feature2 = cuboid2.support_feature(pos21 * -best_sep.1);

    match (&feature1, &feature2) {
        (CuboidFeature::Face(f1), CuboidFeature::Vertex(v2)) => {
            CuboidFeature::face_vertex_contacts(pos12, f1, &best_sep.1, v2, manifold, false)
        }
        (CuboidFeature::Vertex(v1), CuboidFeature::Face(f2)) => {
            CuboidFeature::face_vertex_contacts(&pos12, f2, &-best_sep.1, v1, manifold, true)
        }
        #[cfg(feature = "dim3")]
        (CuboidFeature::Face(f1), CuboidFeature::Edge(e2)) => CuboidFeature::face_edge_contacts(
            &pos12,
            f1,
            &best_sep.1,
            e2,
            prediction,
            manifold,
            false,
        ),
        (CuboidFeature::Face(f1), CuboidFeature::Face(f2)) => CuboidFeature::face_face_contacts(
            &pos12,
            f1,
            &best_sep.1,
            f2,
            prediction,
            manifold,
            false,
        ),
        #[cfg(feature = "dim3")]
        (CuboidFeature::Edge(e1), CuboidFeature::Edge(e2)) => {
            CuboidFeature::edge_edge_contacts(&pos12, e1, &best_sep.1, e2, manifold, false)
        }
        #[cfg(feature = "dim3")]
        (CuboidFeature::Edge(e1), CuboidFeature::Face(f2)) => {
            CuboidFeature::face_edge_contacts(&pos21, f2, &local_n2, e1, prediction, manifold, true)
        }
        _ => unreachable!(), // The other cases are not possible.
    }

    manifold.local_n1 = best_sep.1;
    manifold.local_n2 = local_n2;
    manifold.kinematics.category = KinematicsCategory::PlanePoint;
    manifold.kinematics.radius1 = 0.0;
    manifold.kinematics.radius2 = 0.0;

    // Transfer impulses.
    manifold.match_contacts(&old_manifold_points);
    manifold.sort_contacts(prediction);
}
