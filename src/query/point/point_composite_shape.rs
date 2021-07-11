#![allow(unused_parens)] // Needed by the macro.

use crate::bounding_volume::SimdAABB;
use crate::math::{Point, Real, SimdReal, SIMD_WIDTH};
use crate::partitioning::{SimdBestFirstVisitStatus, SimdBestFirstVisitor};
use crate::query::{
    visitors::CompositePointContainmentTest, PointProjection, PointQuery, PointQueryWithLocation,
};
use crate::shape::{
    Compound, FeatureId, Polyline, SegmentPointLocation, TriMesh, TrianglePointLocation,
    TypedSimdCompositeShape,
};
use na;
use simba::simd::{SimdBool as _, SimdPartialOrd, SimdValue};

impl PointQuery for Polyline {
    #[inline]
    fn project_local_point(&self, point: &Point<Real>, solid: bool) -> PointProjection {
        self.project_local_point_and_get_location(point, solid).0
    }

    #[inline]
    fn project_local_point_and_get_feature(
        &self,
        point: &Point<Real>,
    ) -> (PointProjection, FeatureId) {
        let mut visitor =
            PointCompositeShapeProjWithFeatureBestFirstVisitor::new(self, point, false);
        let (proj, (id, feature)) = self.qbvh().traverse_best_first(&mut visitor).unwrap().1;
        let polyline_feature = self.segment_feature_to_polyline_feature(id, feature);

        (proj, polyline_feature)
    }

    // FIXME: implement distance_to_point too?

    #[inline]
    fn contains_local_point(&self, point: &Point<Real>) -> bool {
        let mut visitor = CompositePointContainmentTest::new(self, point);
        self.qbvh().traverse_depth_first(&mut visitor);
        visitor.found
    }
}

impl PointQuery for TriMesh {
    #[inline]
    fn project_local_point(&self, point: &Point<Real>, solid: bool) -> PointProjection {
        self.project_local_point_and_get_location(point, solid).0
    }

    #[inline]
    fn project_local_point_and_get_feature(
        &self,
        point: &Point<Real>,
    ) -> (PointProjection, FeatureId) {
        let mut visitor =
            PointCompositeShapeProjWithFeatureBestFirstVisitor::new(self, point, false);
        let (proj, (id, _feature)) = self.qbvh().traverse_best_first(&mut visitor).unwrap().1;
        let feature_id = FeatureId::Face(id);
        (proj, feature_id)
    }

    // FIXME: implement distance_to_point too?

    #[inline]
    fn contains_local_point(&self, point: &Point<Real>) -> bool {
        let mut visitor = CompositePointContainmentTest::new(self, point);
        self.qbvh().traverse_depth_first(&mut visitor);
        visitor.found
    }
}

impl PointQuery for Compound {
    #[inline]
    fn project_local_point(&self, point: &Point<Real>, solid: bool) -> PointProjection {
        let mut visitor = PointCompositeShapeProjBestFirstVisitor::new(self, point, solid);
        self.qbvh().traverse_best_first(&mut visitor).unwrap().1 .0
    }

    #[inline]
    fn project_local_point_and_get_feature(
        &self,
        point: &Point<Real>,
    ) -> (PointProjection, FeatureId) {
        (self.project_local_point(point, false), FeatureId::Unknown)
    }

    #[inline]
    fn contains_local_point(&self, point: &Point<Real>) -> bool {
        let mut visitor = CompositePointContainmentTest::new(self, point);
        self.qbvh().traverse_depth_first(&mut visitor);
        visitor.found
    }
}

impl PointQueryWithLocation for Polyline {
    type Location = (u32, SegmentPointLocation);

    #[inline]
    fn project_local_point_and_get_location(
        &self,
        point: &Point<Real>,
        solid: bool,
    ) -> (PointProjection, Self::Location) {
        let mut visitor =
            PointCompositeShapeProjWithLocationBestFirstVisitor::new(self, point, solid);
        self.qbvh().traverse_best_first(&mut visitor).unwrap().1
    }
}

impl PointQueryWithLocation for TriMesh {
    type Location = (u32, TrianglePointLocation);

    #[inline]
    fn project_local_point_and_get_location(
        &self,
        point: &Point<Real>,
        solid: bool,
    ) -> (PointProjection, Self::Location) {
        let mut visitor =
            PointCompositeShapeProjWithLocationBestFirstVisitor::new(self, point, solid);
        self.qbvh().traverse_best_first(&mut visitor).unwrap().1
    }
}

/*
 * Visitors
 */
macro_rules! gen_visitor(
    ($Visitor: ident, $project_local_point: ident, $project_point: ident $(, $Location: ty, $extra_info: ident)* $(| $args: ident)* $(where $PartShapeBound: ident)*) => {
        /// A visitor for the projection of a point on a composite shape.
        pub struct $Visitor<'a, S> {
            shape: &'a S,
            point: &'a Point<Real>,
            simd_point: Point<SimdReal>,
            #[allow(dead_code)] // This won't be used for the projection with feature.
            solid: bool,
        }

        impl<'a, S> $Visitor<'a, S> {
            /// Initialize a visitor for the projection of a point on a composite shape.
            pub fn new(shape: &'a S, point: &'a Point<Real>, solid: bool) -> Self {
                Self {
                    shape,
                    point,
                    simd_point: Point::splat(*point),
                    solid,
                }
            }
        }

        impl<'a, S> SimdBestFirstVisitor<S::PartId, SimdAABB> for $Visitor<'a, S>
        where S: TypedSimdCompositeShape
              $(, $Location: Copy)*
              $(, S::PartShape: $PartShapeBound)* {
            type Result = (PointProjection, (S::PartId $(, $Location)*));

            #[inline]
            fn visit(
                &mut self,
                best: Real,
                aabb: &SimdAABB,
                data: Option<[Option<&S::PartId>; SIMD_WIDTH]>,
            ) -> SimdBestFirstVisitStatus<Self::Result> {
                let dist = aabb.distance_to_local_point(&self.simd_point);
                let mask = dist.simd_lt(SimdReal::splat(best));

                if let Some(data) = data {
                    let mut weights = [0.0; SIMD_WIDTH];
                    let mut results = [None; SIMD_WIDTH];
                    let bitmask = mask.bitmask();

                    for ii in 0..SIMD_WIDTH {
                        if (bitmask & (1 << ii)) != 0 && data[ii].is_some() {
                            let subshape_id = *data[ii].unwrap();
                            self.shape.map_typed_part_at(subshape_id, |part_pos, part_shape| {
                                let (proj $(, $extra_info)*) = if let Some(part_pos) = part_pos {
                                    part_shape.$project_point(
                                        part_pos,
                                        self.point
                                        $(, self.$args)*
                                    )
                                } else {
                                    part_shape.$project_local_point(
                                        self.point
                                        $(, self.$args)*
                                    )
                                };

                                weights[ii] = na::distance(self.point, &proj.point);
                                results[ii] = Some((proj, (subshape_id $(, $extra_info)*)));
                            });
                        }
                    }

                    SimdBestFirstVisitStatus::MaybeContinue {
                        weights: SimdReal::from(weights),
                        mask,
                        results,
                    }
                } else {
                    SimdBestFirstVisitStatus::MaybeContinue {
                        weights: dist,
                        mask,
                        results: [None; SIMD_WIDTH],
                    }
                }
            }
        }
    }
);

gen_visitor!(
    PointCompositeShapeProjBestFirstVisitor,
    project_local_point,
    project_point | solid
);
gen_visitor!(
    PointCompositeShapeProjWithLocationBestFirstVisitor,
    project_local_point_and_get_location,
    project_point_and_get_location,
    <S::PartShape as PointQueryWithLocation>::Location,
    extra_info | solid
    where PointQueryWithLocation
    where Copy
);
gen_visitor!(
    PointCompositeShapeProjWithFeatureBestFirstVisitor,
    project_local_point_and_get_feature,
    project_point_and_get_feature,
    FeatureId,
    extra_info
);
