pub mod query_options_dispatcher;

use core::any::Any;

use crate::math::{Point, Real};
use crate::partitioning::BvhNode;
use crate::query::point::point_composite_shape::query_options_dispatcher::{
    QueryOptionsDispatcher, QueryOptionsDispatcherMap,
};
use crate::query::{
    PointProjection, PointQuery, PointQueryWithLocation, QueryOptions, QueryOptionsNotUsed,
};
use crate::shape::{
    CompositeShapeRef, FeatureId, SegmentPointLocation, TriMesh, TrianglePointLocation,
    TypedCompositeShape,
};
use na;

use crate::shape::{Compound, Polyline};

impl<S: TypedCompositeShape> CompositeShapeRef<'_, S> {
    /// Project a point on this composite shape.
    ///
    /// Returns the projected point as well as the index of the sub-shape of `self` that was hit.
    /// The third tuple element contains some shape-specific information about the projected point.
    #[inline]
    pub fn project_local_point_and_get_location(
        &self,
        point: &Point<Real>,
        max_dist: Real,
        solid: bool,
        options: &dyn QueryOptionsDispatcher,
    ) -> Option<(
        u32,
        (
            PointProjection,
            <S::PartShape as PointQueryWithLocation>::Location,
        ),
    )>
    where
        S::PartShape: PointQueryWithLocation,
    {
        self.0
            .bvh()
            .find_best(
                max_dist,
                |node: &BvhNode, _best_so_far| {
                    node.aabb()
                        .distance_to_local_point(point, true, &QueryOptionsNotUsed)
                },
                |primitive, _best_so_far| {
                    let proj = self.0.map_typed_part_at(primitive, |pose, shape, _| {
                        if let Some(pose) = pose {
                            shape.project_point_and_get_location(
                                pose,
                                point,
                                solid,
                                options.get_option_for_shape(&shape.type_id()),
                            )
                        } else {
                            shape.project_local_point_and_get_location(
                                point,
                                solid,
                                options.get_option_for_shape(&shape.type_id()),
                            )
                        }
                    })?;
                    let cost = na::distance(&proj.0.point, point);
                    Some((cost, proj))
                },
            )
            .map(|(best_id, (_, (proj, location)))| (best_id, (proj, location)))
    }

    /// Project a point on this composite shape.
    ///
    /// Returns the projected point as well as the index of the sub-shape of `self` that was hit.
    /// If `solid` is `false` then the point will be projected to the closest boundary of `self` even
    /// if it is contained by one of its sub-shapes.
    pub fn project_local_point(
        &self,
        point: &Point<Real>,
        solid: bool,
        options: &dyn QueryOptionsDispatcher,
    ) -> (u32, PointProjection) {
        let (best_id, (_, proj)) = self
            .0
            .bvh()
            .find_best(
                Real::MAX,
                |node: &BvhNode, _best_so_far| {
                    node.aabb()
                        .distance_to_local_point(point, true, &QueryOptionsNotUsed)
                },
                |primitive, _best_so_far| {
                    let proj = self.0.map_typed_part_at(primitive, |pose, shape, _| {
                        if let Some(pose) = pose {
                            shape.project_point(
                                pose,
                                point,
                                solid,
                                options.get_option_for_shape(&shape.type_id()),
                            )
                        } else {
                            shape.project_local_point(
                                point,
                                solid,
                                options.get_option_for_shape(&shape.type_id()),
                            )
                        }
                    })?;
                    let dist = na::distance(&proj.point, point);
                    Some((dist, proj))
                },
            )
            .unwrap();
        (best_id, proj)
    }

    /// Project a point on this composite shape.
    ///
    /// Returns the projected point as well as the index of the sub-shape of `self` that was hit.
    /// The third tuple element contains some shape-specific information about the shape feature
    /// hit by the projection.
    #[inline]
    pub fn project_local_point_and_get_feature(
        &self,
        point: &Point<Real>,
        options: &dyn QueryOptionsDispatcher,
    ) -> (u32, (PointProjection, FeatureId)) {
        let (best_id, (_, (proj, feature_id))) = self
            .0
            .bvh()
            .find_best(
                Real::MAX,
                |node: &BvhNode, _best_so_far| {
                    node.aabb()
                        .distance_to_local_point(point, true, &QueryOptionsNotUsed)
                },
                |primitive, _best_so_far| {
                    let proj = self.0.map_typed_part_at(primitive, |pose, shape, _| {
                        if let Some(pose) = pose {
                            shape.project_point_and_get_feature(
                                pose,
                                point,
                                options.get_option_for_shape(&shape.type_id()),
                            )
                        } else {
                            shape.project_local_point_and_get_feature(
                                point,
                                options.get_option_for_shape(&shape.type_id()),
                            )
                        }
                    })?;
                    let cost = na::distance(&proj.0.point, point);
                    Some((cost, proj))
                },
            )
            .unwrap();
        (best_id, (proj, feature_id))
    }

    // TODO: implement distance_to_point too?

    /// Returns the index of any sub-shape of `self` that contains the given point.
    #[inline]
    pub fn contains_local_point(
        &self,
        point: &Point<Real>,
        options: &dyn QueryOptionsDispatcher,
    ) -> Option<u32> {
        self.0
            .bvh()
            .leaves(|node: &BvhNode| node.aabb().contains_local_point(point))
            .find(|leaf_id| {
                self.0
                    .map_typed_part_at(*leaf_id, |pose, shape, _| {
                        if let Some(pose) = pose {
                            shape.contains_point(
                                pose,
                                point,
                                options.get_option_for_shape(&shape.type_id()),
                            )
                        } else {
                            shape.contains_local_point(
                                point,
                                options.get_option_for_shape(&shape.type_id()),
                            )
                        }
                    })
                    .unwrap_or(false)
            })
    }
}

impl PointQuery for Polyline {
    #[inline]
    fn project_local_point(
        &self,
        point: &Point<Real>,
        solid: bool,
        options: &dyn QueryOptions,
    ) -> PointProjection {
        let options = QueryOptionsDispatcherMap::from_dyn_or_default(options);
        CompositeShapeRef(self)
            .project_local_point(point, solid, options)
            .1
    }

    #[inline]
    fn project_local_point_and_get_feature(
        &self,
        point: &Point<Real>,
        options: &dyn QueryOptions,
    ) -> (PointProjection, FeatureId) {
        let options = QueryOptionsDispatcherMap::from_dyn_or_default(options);
        let (seg_id, (proj, feature)) =
            CompositeShapeRef(self).project_local_point_and_get_feature(point, options);
        let polyline_feature = self.segment_feature_to_polyline_feature(seg_id, feature);
        (proj, polyline_feature)
    }

    // TODO: implement distance_to_point too?

    #[inline]
    fn contains_local_point(&self, point: &Point<Real>, options: &dyn QueryOptions) -> bool {
        let options = QueryOptionsDispatcherMap::from_dyn_or_default(options);
        CompositeShapeRef(self)
            .contains_local_point(point, options)
            .is_some()
    }
}

impl PointQuery for TriMesh {
    #[inline]
    fn project_local_point(
        &self,
        point: &Point<Real>,
        solid: bool,
        options: &dyn QueryOptions,
    ) -> PointProjection {
        let options = QueryOptionsDispatcherMap::from_dyn_or_default(options);
        CompositeShapeRef(self)
            .project_local_point(point, solid, options)
            .1
    }

    #[inline]
    fn project_local_point_and_get_feature(
        &self,
        point: &Point<Real>,
        options: &dyn QueryOptions,
    ) -> (PointProjection, FeatureId) {
        #[cfg(feature = "dim3")]
        if self.pseudo_normals().is_some() {
            // If we can, in 3D, take the pseudo-normals into account.
            let (proj, (id, _feature)) =
                self.project_local_point_and_get_location(point, false, options);
            let feature_id = FeatureId::Face(id);
            return (proj, feature_id);
        }
        let options = QueryOptionsDispatcherMap::from_dyn_or_default(options);
        let solid = cfg!(feature = "dim2");
        let (tri_id, proj) = CompositeShapeRef(self).project_local_point(point, solid, options);
        (proj, FeatureId::Face(tri_id))
    }

    // TODO: implement distance_to_point too?

    #[inline]
    fn contains_local_point(&self, point: &Point<Real>, options: &dyn QueryOptions) -> bool {
        #[cfg(feature = "dim3")]
        if self.pseudo_normals.is_some() {
            // If we can, in 3D, take the pseudo-normals into account.
            return self
                .project_local_point_and_get_location(point, true, options)
                .0
                .is_inside;
        }

        let options = QueryOptionsDispatcherMap::from_dyn_or_default(options);
        CompositeShapeRef(self)
            .contains_local_point(point, options)
            .is_some()
    }

    /// Projects a point on `self` transformed by `m`, unless the projection lies further than the given max distance.
    fn project_local_point_with_max_dist(
        &self,
        pt: &Point<Real>,
        solid: bool,
        max_dist: Real,
        options: &dyn QueryOptions,
    ) -> Option<PointProjection> {
        self.project_local_point_and_get_location_with_max_dist(pt, solid, max_dist, options)
            .map(|proj| proj.0)
    }
}

impl PointQuery for Compound {
    #[inline]
    fn project_local_point(
        &self,
        point: &Point<Real>,
        solid: bool,
        options: &dyn QueryOptions,
    ) -> PointProjection {
        let options = QueryOptionsDispatcherMap::from_dyn_or_default(options);
        CompositeShapeRef(self)
            .project_local_point(point, solid, options)
            .1
    }

    #[inline]
    fn project_local_point_and_get_feature(
        &self,
        point: &Point<Real>,
        options: &dyn QueryOptions,
    ) -> (PointProjection, FeatureId) {
        let options = QueryOptionsDispatcherMap::from_dyn_or_default(options);
        (
            CompositeShapeRef(self)
                .project_local_point_and_get_feature(point, options)
                .1
                 .0,
            FeatureId::Unknown,
        )
    }

    #[inline]
    fn contains_local_point(&self, point: &Point<Real>, options: &dyn QueryOptions) -> bool {
        let options = QueryOptionsDispatcherMap::from_dyn_or_default(options);
        CompositeShapeRef(self)
            .contains_local_point(point, options)
            .is_some()
    }
}

impl PointQueryWithLocation for Polyline {
    type Location = (u32, SegmentPointLocation);

    #[inline]
    fn project_local_point_and_get_location(
        &self,
        point: &Point<Real>,
        solid: bool,
        options: &dyn QueryOptions,
    ) -> (PointProjection, Self::Location) {
        let options = QueryOptionsDispatcherMap::from_dyn_or_default(options);
        let (seg_id, (proj, loc)) = CompositeShapeRef(self)
            .project_local_point_and_get_location(point, Real::MAX, solid, options)
            .unwrap();
        (proj, (seg_id, loc))
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
        options: &dyn QueryOptions,
    ) -> (PointProjection, Self::Location) {
        self.project_local_point_and_get_location_with_max_dist(point, solid, Real::MAX, options)
            .unwrap()
    }

    /// Projects a point on `self`, with a maximum projection distance.
    fn project_local_point_and_get_location_with_max_dist(
        &self,
        point: &Point<Real>,
        solid: bool,
        max_dist: Real,
        options: &dyn QueryOptions,
    ) -> Option<(PointProjection, Self::Location)> {
        let options = QueryOptionsDispatcherMap::from_dyn_or_default(options);
        #[allow(unused_mut)] // mut is needed in 3D.
        if let Some((part_id, (mut proj, location))) = CompositeShapeRef(self)
            .project_local_point_and_get_location(point, max_dist, solid, options)
        {
            #[cfg(feature = "dim3")]
            if let Some(pseudo_normals) = self.pseudo_normals_if_oriented() {
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
