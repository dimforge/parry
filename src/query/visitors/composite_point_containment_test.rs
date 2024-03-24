use crate::bounding_volume::SimdAabb;
use crate::math::{Point, Real, SimdReal, SIMD_WIDTH};
use crate::partitioning::{SimdVisitStatus, SimdVisitor};
use crate::query::point::point_query::PointQuery;
use crate::shape::TypedSimdCompositeShape;
use crate::utils::IsometryOpt;
use simba::simd::{SimdBool as _, SimdValue};

/// Visitor for checking if a composite shape contains a specific point.
pub struct CompositePointContainmentTest<'a, S: 'a> {
    /// The composite shape on which the point containment test should be performed.
    pub shape: &'a S,
    /// The point to be tested.
    pub point: &'a Point<Real>,
    /// A traversal will set this to `true` if the point is inside of `self.shape`.
    pub found: bool,
}

impl<'a, S> CompositePointContainmentTest<'a, S> {
    /// Creates a new visitor for the testing containment of the given `point`
    /// into the given `shape`.
    pub fn new(shape: &'a S, point: &'a Point<Real>) -> Self {
        Self {
            shape,
            point,
            found: false,
        }
    }
}

impl<'a, S: TypedSimdCompositeShape> SimdVisitor<S::PartId, SimdAabb>
    for CompositePointContainmentTest<'a, S>
{
    #[inline]
    fn visit(
        &mut self,
        bv: &SimdAabb,
        b: Option<[Option<&S::PartId>; SIMD_WIDTH]>,
    ) -> SimdVisitStatus {
        let simd_point: Point<SimdReal> = Point::splat(*self.point);
        let mask = bv.contains_local_point(&simd_point);

        if let Some(data) = b {
            let bitmask = mask.bitmask();

            for (ii, data) in data.into_iter().enumerate() {
                if (bitmask & (1 << ii)) != 0 {
                    let Some(data) = data else { continue };
                    self.shape.map_typed_part_at(*data, |part_pos, obj| {
                        if obj.contains_local_point(&part_pos.inverse_transform_point(self.point)) {
                            self.found = true;
                        }
                    });

                    if self.found {
                        return SimdVisitStatus::ExitEarly;
                    }
                }
            }
        }

        SimdVisitStatus::MaybeContinue(mask)
    }
}
