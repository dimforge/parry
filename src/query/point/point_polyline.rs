use crate::bounding_volume::SimdAABB;
use crate::math::{Point, Real, SimdReal, SIMD_WIDTH};
use crate::partitioning::{SimdBestFirstVisitStatus, SimdBestFirstVisitor};
use crate::query::{
    visitors::CompositePointContainmentTest, PointProjection, PointQuery, PointQueryWithLocation,
};
use crate::shape::{FeatureId, Polyline, SegmentPointLocation};
use na;
use simba::simd::{SimdBool as _, SimdPartialOrd, SimdValue};

impl PointQuery for Polyline {
    #[inline]
    fn project_local_point(&self, point: &Point<Real>, solid: bool) -> PointProjection {
        self.project_local_point_with_location(point, solid).0
    }

    #[inline]
    fn project_local_point_and_get_feature(
        &self,
        pt: &Point<Real>,
    ) -> (PointProjection, FeatureId) {
        let mut visitor = PolylinePointProjWithFeatureVisitor {
            polyline: self,
            point: pt,
            simd_point: Point::splat(*pt),
        };

        let (proj, (id, feature)) = self.quadtree().traverse_best_first(&mut visitor).unwrap().1;
        let polyline_feature = self.segment_feature_to_polyline_feature(id, feature);

        (proj, polyline_feature)
    }

    // FIXME: implement distance_to_point too?

    #[inline]
    fn contains_local_point(&self, point: &Point<Real>) -> bool {
        let mut visitor = CompositePointContainmentTest {
            shape: self,
            point,
            found: false,
        };

        self.quadtree().traverse_depth_first(&mut visitor);

        visitor.found
    }
}

impl PointQueryWithLocation for Polyline {
    type Location = (u32, SegmentPointLocation);

    #[inline]
    fn project_local_point_with_location(
        &self,
        point: &Point<Real>,
        _: bool,
    ) -> (PointProjection, Self::Location) {
        let mut visitor = PolylinePointProjVisitor {
            polyline: self,
            point,
            simd_point: Point::splat(*point),
        };

        self.quadtree().traverse_best_first(&mut visitor).unwrap().1
    }
}

/*
 * Visitors
 */
macro_rules! gen_visitor(
    ($Visitor: ident, $Location: ty, $project: ident $(, $args: ident)*) => {
        struct $Visitor<'a> {
            polyline: &'a Polyline,
            point: &'a Point<Real>,
            simd_point: Point<SimdReal>,
        }

        impl<'a> SimdBestFirstVisitor<u32, SimdAABB> for $Visitor<'a> {
            type Result = (PointProjection, (u32, $Location));

            #[inline]
            fn visit(
                &mut self,
                best: Real,
                aabb: &SimdAABB,
                data: Option<[Option<&u32>; SIMD_WIDTH]>,
            ) -> SimdBestFirstVisitStatus<Self::Result> {
                let dist = aabb.distance_to_local_point(&self.simd_point);
                let mask = dist.simd_lt(SimdReal::splat(best));

                if let Some(data) = data {
                    let mut weights = [0.0; SIMD_WIDTH];
                    let mut results = [None; SIMD_WIDTH];
                    let bitmask = mask.bitmask();

                    for ii in 0..SIMD_WIDTH {
                        if (bitmask & (1 << ii)) != 0 && data[ii].is_some(){
                            let (proj, extra_info) = self.polyline.segment(*data[ii].unwrap()).$project(
                                self.point
                                $(, $args)*
                            );

                            let extra_info = (*data[ii].unwrap(), extra_info);
                            weights[ii] = na::distance(self.point, &proj.point);
                            results[ii] = Some((proj, extra_info));
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
    PolylinePointProjVisitor,
    SegmentPointLocation,
    project_local_point_with_location,
    true
);
gen_visitor!(
    PolylinePointProjWithFeatureVisitor,
    FeatureId,
    project_local_point_and_get_feature
);
