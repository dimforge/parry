#[cfg(feature = "dim2")]
use crate::math::Vector;
use crate::math::{Isometry, Real};
use crate::query::{sat, ContactManifold};
use crate::shape::PolygonalFeature;
use crate::shape::{Cuboid, Shape, Triangle};

/// Computes the contact manifold between a cuboid and a triangle represented as `Shape` trait-objects.
pub fn contact_manifold_cuboid_triangle_shapes<ManifoldData, ContactData>(
    pos12: &Isometry<Real>,
    shape1: &dyn Shape,
    shape2: &dyn Shape,
    prediction: Real,
    manifold: &mut ContactManifold<ManifoldData, ContactData>,
) where
    ContactData: Default + Copy,
{
    if let (Some(cuboid1), Some(triangle2)) = (shape1.as_cuboid(), shape2.as_triangle()) {
        contact_manifold_cuboid_triangle(
            pos12,
            &pos12.inverse(),
            cuboid1,
            triangle2,
            prediction,
            manifold,
            false,
        );
    } else if let (Some(triangle1), Some(cuboid2)) = (shape1.as_triangle(), shape2.as_cuboid()) {
        contact_manifold_cuboid_triangle(
            &pos12.inverse(),
            pos12,
            cuboid2,
            triangle1,
            prediction,
            manifold,
            true,
        );
    }
}

/// Computes the contact manifold between a cuboid and a triangle.
pub fn contact_manifold_cuboid_triangle<'a, ManifoldData, ContactData>(
    pos12: &Isometry<Real>,
    pos21: &Isometry<Real>,
    cuboid1: &'a Cuboid,
    triangle2: &'a Triangle,
    prediction: Real,
    manifold: &mut ContactManifold<ManifoldData, ContactData>,
    flipped: bool,
) where
    ContactData: Default + Copy,
{
    if (!flipped && manifold.try_update_contacts(&pos12))
        || (flipped && manifold.try_update_contacts(&pos21))
    {
        return;
    }

    /*
     *
     * Point-Face cases.
     *
     */
    let sep1 =
        sat::cuboid_support_map_find_local_separating_normal_oneway(cuboid1, triangle2, &pos12);
    if sep1.0 > prediction {
        manifold.clear();
        return;
    }

    let sep2 = sat::triangle_cuboid_find_local_separating_normal_oneway(triangle2, cuboid1, &pos21);
    if sep2.0 > prediction {
        manifold.clear();
        return;
    }

    /*
     *
     * Edge-Edge cases.
     *
     */
    #[cfg(feature = "dim2")]
    let sep3 = (-Real::MAX, Vector::x()); // This case does not exist in 2D.
    #[cfg(feature = "dim3")]
    let sep3 = sat::cuboid_triangle_find_local_separating_edge_twoway(cuboid1, triangle2, &pos12);
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

    let feature1;
    let feature2;
    let normal2 = pos21 * -best_sep.1;

    #[cfg(feature = "dim2")]
    {
        feature1 = cuboid1.support_face(best_sep.1);
        feature2 = triangle2.support_face(normal2);
    }
    #[cfg(feature = "dim3")]
    {
        feature1 = cuboid1.support_face(best_sep.1);
        feature2 = PolygonalFeature::from(*triangle2);
    }

    // We do this clone to perform contact tracking and transfer impulses.
    // FIXME: find a more efficient way of doing this.
    let old_manifold_points = manifold.points.clone();
    manifold.clear();

    PolygonalFeature::contacts(
        pos12,
        pos21,
        &best_sep.1,
        &normal2,
        &feature1,
        &feature2,
        prediction,
        manifold,
        flipped,
    );

    if flipped {
        manifold.local_n1 = normal2;
        manifold.local_n2 = best_sep.1;
    } else {
        manifold.local_n1 = best_sep.1;
        manifold.local_n2 = normal2;
    }

    // Transfer impulses.
    manifold.match_contacts(&old_manifold_points);
}
