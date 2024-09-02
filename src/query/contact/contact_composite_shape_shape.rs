use crate::bounding_volume::BoundingVolume;
use crate::math::{Isometry, Real};
use crate::query::visitors::BoundingVolumeIntersectionsVisitor;
use crate::query::{Contact, QueryDispatcher};
use crate::shape::{Shape, SimdCompositeShape};
use crate::utils::IsometryOpt;

/// Best contact between a composite shape (`Mesh`, `Compound`) and any other shape.
pub fn contact_composite_shape_shape<D, G1>(
    dispatcher: &D,
    pos12: &Isometry<Real>,
    g1: &G1,
    g2: &dyn Shape,
    prediction: Real,
) -> Option<Contact>
where
    D: ?Sized + QueryDispatcher,
    G1: ?Sized + SimdCompositeShape,
{
    // Find new collisions
    let ls_aabb2 = g2.compute_aabb(pos12).loosened(prediction);
    let mut res = None::<Contact>;

    let mut leaf_callback = |i: &_| {
        g1.map_part_at(*i, &mut |part_pos1, part1, _| {
            if let Ok(Some(mut c)) =
                dispatcher.contact(&part_pos1.inv_mul(pos12), part1, g2, prediction)
            {
                let replace = res.map_or(true, |cbest| c.dist < cbest.dist);

                if replace {
                    if let Some(part_pos1) = part_pos1 {
                        c.transform1_by_mut(part_pos1);
                    }
                    res = Some(c)
                }
            }
        });

        true
    };

    let mut visitor = BoundingVolumeIntersectionsVisitor::new(&ls_aabb2, &mut leaf_callback);
    let _ = g1.qbvh().traverse_depth_first(&mut visitor);
    res
}

/// Best contact between a shape and a composite (`Mesh`, `Compound`) shape.
pub fn contact_shape_composite_shape<D, G2>(
    dispatcher: &D,
    pos12: &Isometry<Real>,
    g1: &dyn Shape,
    g2: &G2,
    prediction: Real,
) -> Option<Contact>
where
    D: ?Sized + QueryDispatcher,
    G2: ?Sized + SimdCompositeShape,
{
    contact_composite_shape_shape(dispatcher, &pos12.inverse(), g2, g1, prediction)
        .map(|c| c.flipped())
}
