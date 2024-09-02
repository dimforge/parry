#![allow(unused_parens)] // Needed by the macro.

use crate::bounding_volume::SimdAabb;
use crate::math::{Point, Real, SimdReal, SIMD_WIDTH};
use crate::partitioning::{SimdBestFirstVisitStatus, SimdBestFirstVisitor};
use crate::query::visitors::CompositePointContainmentTest;
use crate::query::{PointProjection, PointQuery, PointQueryWithLocation};
use crate::shape::{
    FeatureId, SegmentPointLocation, TriMesh, TrianglePointLocation, TypedSimdCompositeShape,
};
use na;
use simba::simd::{SimdBool as _, SimdPartialOrd, SimdValue};

use crate::shape::{Compound, Polyline};

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

    // TODO: implement distance_to_point too?

    #[inline]
    fn contains_local_point(&self, point: &Point<Real>) -> bool {
        let mut visitor = CompositePointContainmentTest::new(self, point);
        let _ = self.qbvh().traverse_depth_first(&mut visitor);
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
        #[cfg(feature = "dim3")]
        if self.pseudo_normals().is_some() {
            // If we can, in 3D, take the pseudo-normals into account.
            let (proj, (id, _feature)) = self.project_local_point_and_get_location(point, false);
            let feature_id = FeatureId::Face(id);
            return (proj, feature_id);
        }

        let solid = cfg!(feature = "dim2");

        let mut visitor =
            PointCompositeShapeProjWithFeatureBestFirstVisitor::new(self, point, solid);
        let (proj, (id, _feature)) = self.qbvh().traverse_best_first(&mut visitor).unwrap().1;
        let feature_id = FeatureId::Face(id);
        (proj, feature_id)
    }

    // TODO: implement distance_to_point too?

    #[inline]
    fn contains_local_point(&self, point: &Point<Real>) -> bool {
        #[cfg(feature = "dim3")]
        if self.pseudo_normals.is_some() {
            // If we can, in 3D, take the pseudo-normals into account.
            return self
                .project_local_point_and_get_location(point, true)
                .0
                .is_inside;
        }

        let mut visitor = CompositePointContainmentTest::new(self, point);
        let _ = self.qbvh().traverse_depth_first(&mut visitor);
        visitor.found
    }

    /// Projects a point on `self` transformed by `m`, unless the projection lies further than the given max distance.
    fn project_local_point_with_max_dist(
        &self,
        pt: &Point<Real>,
        solid: bool,
        max_dist: Real,
    ) -> Option<PointProjection> {
        self.project_local_point_and_get_location_with_max_dist(pt, solid, max_dist)
            .map(|proj| proj.0)
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
        let _ = self.qbvh().traverse_depth_first(&mut visitor);
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
    #[allow(unused_mut)] // Because we need mut in 3D but not in 2D.
    fn project_local_point_and_get_location(
        &self,
        point: &Point<Real>,
        solid: bool,
    ) -> (PointProjection, Self::Location) {
        self.project_local_point_and_get_location_with_max_dist(point, solid, Real::MAX)
            .unwrap()
    }

    /// Projects a point on `self`, with a maximum projection distance.
    fn project_local_point_and_get_location_with_max_dist(
        &self,
        point: &Point<Real>,
        solid: bool,
        max_dist: Real,
    ) -> Option<(PointProjection, Self::Location)> {
        let mut visitor =
            PointCompositeShapeProjWithLocationBestFirstVisitor::new(self, point, solid);

        #[allow(unused_mut)] // mut is needed in 3D.
        if let Some((_, (mut proj, (part_id, location)))) =
            self.qbvh()
                .traverse_best_first_node(&mut visitor, 0, max_dist)
        {
            #[cfg(feature = "dim3")]
            if let Some(pseudo_normals) = self.pseudo_normals() {
                let pseudo_normal = match location {
                    TrianglePointLocation::OnFace(..) | TrianglePointLocation::OnSolid => {
                        Some(self.triangle(part_id).scaled_normal())
                    }
                    TrianglePointLocation::OnEdge(i, _) => pseudo_normals
                        .edges_pseudo_normal
                        .get(part_id as usize)
                        .map(|pn| pn[i as usize]),
                    TrianglePointLocation::OnVertex(i) => {
                        let idx = self.indices()[part_id as usize];
                        pseudo_normals
                            .vertices_pseudo_normal
                            .get(idx[i as usize] as usize)
                            .copied()
                    }
                };

                if let Some(pseudo_normal) = pseudo_normal {
                    let dpt = point - proj.point;
                    proj.is_inside = dpt.dot(&pseudo_normal) <= 0.0;
                }
            }

            Some((proj, (part_id, location)))
        } else {
            None
        }
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

        impl<'a, S> SimdBestFirstVisitor<S::PartId, SimdAabb> for $Visitor<'a, S>
        where S: TypedSimdCompositeShape
              $(, $Location: Copy)*
              $(, S::PartShape: $PartShapeBound)* {
            type Result = (PointProjection, (S::PartId $(, $Location)*));

            #[inline]
            fn visit(
                &mut self,
                best: Real,
                aabb: &SimdAabb,
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
                            let mut is_inside = false;
                            let subshape_id = *data[ii].unwrap();
                            self.shape.map_typed_part_at(subshape_id, |part_pos, part_shape, _| {
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

                                is_inside = proj.is_inside;
                                weights[ii] = na::distance(self.point, &proj.point);
                                results[ii] = Some((proj, (subshape_id $(, $extra_info)*)));
                            });

                            if self.solid && is_inside {
                                return SimdBestFirstVisitStatus::ExitEarly(results[ii]);
                            }
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
