use crate::bounding_volume::SimdAabb;
use crate::math::{Isometry, Real, SIMD_WIDTH};
use crate::partitioning::{SimdVisitStatus, SimdVisitor};
use crate::query::QueryDispatcher;
use crate::shape::{Shape, TypedSimdCompositeShape};
use crate::utils::IsometryOpt;
use simba::simd::SimdBool as _;

/// Intersection test between a composite shape (`Mesh`, `Compound`) and any other shape.
pub fn intersection_test_composite_shape_shape<D, G1>(
    dispatcher: &D,
    pos12: &Isometry<Real>,
    g1: &G1,
    g2: &dyn Shape,
) -> bool
where
    D: ?Sized + QueryDispatcher,
    G1: ?Sized + TypedSimdCompositeShape,
{
    let mut visitor = IntersectionCompositeShapeShapeVisitor::new(dispatcher, pos12, g1, g2);

    let _ = g1.typed_qbvh().traverse_depth_first(&mut visitor);
    visitor.found_intersection
}

/// Proximity between a shape and a composite (`Mesh`, `Compound`) shape.
pub fn intersection_test_shape_composite_shape<D, G2>(
    dispatcher: &D,
    pos12: &Isometry<Real>,
    g1: &dyn Shape,
    g2: &G2,
) -> bool
where
    D: ?Sized + QueryDispatcher,
    G2: ?Sized + TypedSimdCompositeShape,
{
    intersection_test_composite_shape_shape(dispatcher, &pos12.inverse(), g2, g1)
}

/// A visitor for checking if a composite-shape and a shape intersect.
pub struct IntersectionCompositeShapeShapeVisitor<'a, D: ?Sized, G1: ?Sized + 'a> {
    ls_aabb2: SimdAabb,

    dispatcher: &'a D,
    pos12: &'a Isometry<Real>,
    g1: &'a G1,
    g2: &'a dyn Shape,

    found_intersection: bool,
}

impl<'a, D, G1> IntersectionCompositeShapeShapeVisitor<'a, D, G1>
where
    D: ?Sized + QueryDispatcher,
    G1: ?Sized + TypedSimdCompositeShape,
{
    /// Initialize a visitor for checking if a composite-shape and a shape intersect.
    pub fn new(
        dispatcher: &'a D,
        pos12: &'a Isometry<Real>,
        g1: &'a G1,
        g2: &'a dyn Shape,
    ) -> IntersectionCompositeShapeShapeVisitor<'a, D, G1> {
        let ls_aabb2 = g2.compute_aabb(pos12);

        IntersectionCompositeShapeShapeVisitor {
            dispatcher,
            ls_aabb2: SimdAabb::splat(ls_aabb2),
            pos12,
            g1,
            g2,
            found_intersection: false,
        }
    }
}

impl<'a, D, G1> SimdVisitor<G1::PartId, SimdAabb>
    for IntersectionCompositeShapeShapeVisitor<'a, D, G1>
where
    D: ?Sized + QueryDispatcher,
    G1: ?Sized + TypedSimdCompositeShape,
{
    fn visit(
        &mut self,
        bv: &SimdAabb,
        data: Option<[Option<&G1::PartId>; SIMD_WIDTH]>,
    ) -> SimdVisitStatus {
        let mask = self.ls_aabb2.intersects(bv);

        if let Some(data) = data {
            let bitmask = mask.bitmask();
            let mut found_intersection = false;

            for (ii, data) in data.into_iter().enumerate() {
                if (bitmask & (1 << ii)) != 0 {
                    let Some(data) = data else { continue };
                    let part_id = *data;
                    self.g1.map_untyped_part_at(part_id, |part_pos1, g1, _| {
                        found_intersection = self.dispatcher.intersection_test(
                            &part_pos1.inv_mul(self.pos12),
                            g1,
                            self.g2,
                        ) == Ok(true);
                    });

                    if found_intersection {
                        self.found_intersection = true;
                        return SimdVisitStatus::ExitEarly;
                    }
                }
            }
        }

        SimdVisitStatus::MaybeContinue(mask)
    }
}
