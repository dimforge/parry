use crate::bounding_volume::BoundingVolume;
use crate::math::{Isometry, Real};
use crate::query::visitors::BoundingVolumeIntersectionsVisitor;
use crate::query::{Contact, QueryDispatcher};
use crate::shape::{Shape, SimdCompositeShape};
use crate::utils::IsometryOpt;

/// Best contact between a composite shape (`Mesh`, `Compound`) and any other shape.
pub fn contact_composite_shape_shape<D: ?Sized, G1: ?Sized>(
    dispatcher: &D,
    pos12: &Isometry<Real>,
    g1: &G1,
    g2: &dyn Shape,
    prediction: Real,
) -> Option<Contact>
where
    D: QueryDispatcher,
    G1: SimdCompositeShape,
{
    // Find new collisions
    let ls_aabb2 = g2.compute_aabb(pos12).loosened(prediction);
    let mut res = None::<Contact>;

    let mut leaf_callback = |i: &_| {
        g1.map_part_at(*i, &mut |part_pos1, part1| {
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
    g1.qbvh().traverse_depth_first(&mut visitor);
    res
}

/// Best contact between a shape and a composite (`Mesh`, `Compound`) shape.
pub fn contact_shape_composite_shape<D: ?Sized, G2: ?Sized>(
    dispatcher: &D,
    pos12: &Isometry<Real>,
    g1: &dyn Shape,
    g2: &G2,
    prediction: Real,
) -> Option<Contact>
where
    D: QueryDispatcher,
    G2: SimdCompositeShape,
{
    contact_composite_shape_shape(dispatcher, &pos12.inverse(), g2, g1, prediction)
        .map(|c| c.flipped())
}
