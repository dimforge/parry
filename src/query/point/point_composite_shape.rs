#![allow(unused_parens)] // Needed by the macro.

use crate::math::{Real, Vector};
use crate::partitioning::BvhNode;
use crate::query::{PointProjection, PointQuery, PointQueryWithLocation};
use crate::shape::{
    CompositeShapeRef, FeatureId, SegmentPointLocation, SubShapeId, TriMesh, TrianglePointLocation,
    TypedCompositeShape,
};

use crate::shape::{Compound, Polyline};

/// The feature of a triangle a point projected onto, as `Triangle`'s own point query reports it.
fn triangle_point_location_feature(location: TrianglePointLocation) -> FeatureId {
    match location {
        TrianglePointLocation::OnVertex(i) => FeatureId::Vertex(i),
        #[cfg(feature = "dim3")]
        TrianglePointLocation::OnEdge(i, _) => FeatureId::Edge(i),
        #[cfg(feature = "dim2")]
        TrianglePointLocation::OnEdge(i, _) => FeatureId::Face(i),
        TrianglePointLocation::OnFace(i, _) => FeatureId::Face(i),
        TrianglePointLocation::OnSolid => FeatureId::Face(0),
    }
}

impl<S: TypedCompositeShape> CompositeShapeRef<'_, S> {
    /// Project a point on this composite shape.
    ///
    /// The projection's `subshape` says which sub-shape of `self` answered. The second tuple
    /// element contains some shape-specific information about the projected point.
    #[inline]
    pub fn project_local_point_and_get_location(
        &self,
        point: Vector,
        max_dist: Real,
        solid: bool,
    ) -> Option<(
        PointProjection,
        <S::PartShape as PointQueryWithLocation>::Location,
    )>
    where
        S::PartShape: PointQueryWithLocation,
    {
        self.0
            .bvh()
            .find_best(
                max_dist,
                |node: &BvhNode, _best_so_far| node.aabb().distance_to_local_point(point, true),
                |primitive, _best_so_far| {
                    let proj = self.0.map_typed_part_at(primitive, |pose, shape, _| {
                        if let Some(pose) = pose {
                            shape.project_point_and_get_location(pose, point, solid)
                        } else {
                            shape.project_local_point_and_get_location(point, solid)
                        }
                    })?;
                    let cost = (proj.0.point - point).length();
                    Some((cost, proj))
                },
            )
            .map(|(best_id, (_, (proj, location)))| (proj.with_subshape(best_id), location))
    }

    /// Project a point on this composite shape.
    ///
    /// The projection's `subshape` says which sub-shape of `self` answered. If `solid` is `false`
    /// then the point will be projected to the closest boundary of `self` even if it is contained
    /// by one of its sub-shapes.
    pub fn project_local_point(
        &self,
        point: Vector,
        max_dist: Real,
        solid: bool,
    ) -> Option<PointProjection> {
        let (best_id, (_, proj)) = self.0.bvh().find_best(
            max_dist,
            |node: &BvhNode, _best_so_far| node.aabb().distance_to_local_point(point, true),
            |primitive, _best_so_far| {
                let proj = self.0.map_typed_part_at(primitive, |pose, shape, _| {
                    if let Some(pose) = pose {
                        shape.project_point(pose, point, solid)
                    } else {
                        shape.project_local_point(point, solid)
                    }
                })?;
                let dist = (proj.point - point).length();
                Some((dist, proj))
            },
        )?;
        Some(proj.with_subshape(best_id))
    }

    /// Project a point on this composite shape.
    ///
    /// The projection's `subshape` says which sub-shape of `self` answered. The second tuple
    /// element is the feature of that sub-shape the projection landed on.
    #[inline]
    pub fn project_local_point_and_get_feature(
        &self,
        point: Vector,
        max_dist: Real,
    ) -> Option<(PointProjection, FeatureId)> {
        let (best_id, (_, (proj, feature_id))) = self.0.bvh().find_best(
            max_dist,
            |node: &BvhNode, _best_so_far| node.aabb().distance_to_local_point(point, true),
            |primitive, _best_so_far| {
                let proj = self.0.map_typed_part_at(primitive, |pose, shape, _| {
                    if let Some(pose) = pose {
                        shape.project_point_and_get_feature(pose, point)
                    } else {
                        shape.project_local_point_and_get_feature(point)
                    }
                })?;
                let cost = (proj.0.point - point).length();
                Some((cost, proj))
            },
        )?;
        Some((proj.with_subshape(best_id), feature_id))
    }

    // TODO: implement distance_to_point too?

    /// Returns the index of any sub-shape of `self` that contains the given point.
    #[inline]
    pub fn contains_local_point(&self, point: Vector) -> Option<SubShapeId> {
        self.0
            .bvh()
            .leaves(|node: &BvhNode| node.aabb().contains_local_point(point))
            .find(|leaf_id| {
                self.0
                    .map_typed_part_at(*leaf_id, |pose, shape, _| {
                        if let Some(pose) = pose {
                            shape.contains_point(pose, point)
                        } else {
                            shape.contains_local_point(point)
                        }
                    })
                    .unwrap_or(false)
            })
    }
}

impl PointQuery for Polyline {
    #[inline]
    fn project_local_point(&self, point: Vector, solid: bool) -> PointProjection {
        self.project_local_point_and_get_location(point, solid).0
    }

    #[inline]
    #[allow(unused_mut)] // Because we need mut in 2D but not in 3D.
    fn project_local_point_and_get_feature(&self, point: Vector) -> (PointProjection, FeatureId) {
        // Every comparison involving a NaN is false, so the traversal finds no candidate
        // at all when `point` (or `self`) isn’t finite. Report `point` itself rather than
        // an arbitrary projection onto whichever part we happened to pick.
        let Some((mut proj, feature)) =
            CompositeShapeRef(self).project_local_point_and_get_feature(point, Real::MAX)
        else {
            return (PointProjection::new(false, point), FeatureId::Unknown);
        };

        // A point behind the outward pseudo-normal is inside.
        #[cfg(feature = "dim2")]
        if let Some(constraints) = self.segment_normal_constraints(proj.subshape) {
            let pseudo_normal = match feature {
                FeatureId::Vertex(i) => constraints.edges[i as usize],
                _ => constraints.face,
            };
            proj.is_inside = (point - proj.point).dot(pseudo_normal) <= 0.0;
        }

        // The feature is the segment's own; `proj.subshape` says which segment it belongs to.
        (proj, feature)
    }

    // TODO: implement distance_to_point too?

    #[inline]
    fn contains_local_point(&self, point: Vector) -> bool {
        // An oriented polyline has a solid interior; reuse the projection's inside test.
        #[cfg(feature = "dim2")]
        if self.flags().contains(crate::shape::PolylineFlags::ORIENTED) {
            return self
                .project_local_point_and_get_location(point, true)
                .0
                .is_inside;
        }

        CompositeShapeRef(self)
            .contains_local_point(point)
            .is_some()
    }
}

impl PointQuery for TriMesh {
    #[inline]
    fn project_local_point(&self, point: Vector, solid: bool) -> PointProjection {
        self.project_local_point_with_max_dist(point, solid, Real::MAX)
            // Shouldn’t happen (trimesh must not be empty). But return something
            // instead of crashing with `unwrap`.
            .unwrap_or(PointProjection::new(false, point))
    }

    #[inline]
    fn project_local_point_and_get_feature(&self, point: Vector) -> (PointProjection, FeatureId) {
        #[cfg(feature = "dim3")]
        if self.pseudo_normals().is_some() {
            // If we can, in 3D, take the pseudo-normals into account. The location carries the
            // triangle's own feature; `proj.subshape` says which triangle it belongs to.
            let (proj, (_, location)) = self.project_local_point_and_get_location(point, false);
            return (proj, triangle_point_location_feature(location));
        }

        let solid = cfg!(feature = "dim2");
        // No candidate: `point` (or `self`) isn’t finite. See
        // `Polyline::project_local_point_and_get_feature`.
        let Some((proj, location)) =
            CompositeShapeRef(self).project_local_point_and_get_location(point, Real::MAX, solid)
        else {
            return (PointProjection::new(false, point), FeatureId::Unknown);
        };
        // The feature is the triangle's own; `proj.subshape` says which triangle it belongs to.
        (proj, triangle_point_location_feature(location))
    }

    // TODO: implement distance_to_point too?

    #[inline]
    fn contains_local_point(&self, point: Vector) -> bool {
        #[cfg(feature = "dim3")]
        if self.pseudo_normals.is_some() {
            // If we can, in 3D, take the pseudo-normals into account.
            return self
                .project_local_point_and_get_location(point, true)
                .0
                .is_inside;
        }

        CompositeShapeRef(self)
            .contains_local_point(point)
            .is_some()
    }

    /// Projects a point on `self` transformed by `m`, unless the projection lies further than the given max distance.
    fn project_local_point_with_max_dist(
        &self,
        pt: Vector,
        solid: bool,
        max_dist: Real,
    ) -> Option<PointProjection> {
        self.project_local_point_and_get_location_with_max_dist(pt, solid, max_dist)
            .map(|proj| proj.0)
    }
}

impl PointQuery for Compound {
    #[inline]
    fn project_local_point(&self, point: Vector, solid: bool) -> PointProjection {
        CompositeShapeRef(self)
            .project_local_point(point, Real::MAX, solid)
            // No candidate: `point` (or `self`) isn’t finite. See
            // `Polyline::project_local_point_and_get_feature`.
            .unwrap_or(PointProjection::new(false, point))
    }

    #[inline]
    fn project_local_point_and_get_feature(&self, point: Vector) -> (PointProjection, FeatureId) {
        // The feature is the part's own; `proj.subshape` says which part it belongs to.
        CompositeShapeRef(self)
            .project_local_point_and_get_feature(point, Real::MAX)
            // No candidate: `point` (or `self`) isn’t finite. See
            // `Polyline::project_local_point_and_get_feature`.
            .unwrap_or((PointProjection::new(false, point), FeatureId::Unknown))
    }

    #[inline]
    fn contains_local_point(&self, point: Vector) -> bool {
        CompositeShapeRef(self)
            .contains_local_point(point)
            .is_some()
    }
}

impl PointQueryWithLocation for Polyline {
    type Location = (u32, SegmentPointLocation);

    #[inline]
    fn project_local_point_and_get_location(
        &self,
        point: Vector,
        solid: bool,
    ) -> (PointProjection, Self::Location) {
        self.project_local_point_and_get_location_with_max_dist(point, solid, Real::MAX)
            // No candidate: `point` (or `self`) isn’t finite. See
            // `Polyline::project_local_point_and_get_feature`.
            .unwrap_or((
                PointProjection::new(false, point),
                (0, SegmentPointLocation::OnVertex(0)),
            ))
    }

    /// Projects a point on `self`, with a maximum projection distance.
    fn project_local_point_and_get_location_with_max_dist(
        &self,
        point: Vector,
        solid: bool,
        max_dist: Real,
    ) -> Option<(PointProjection, Self::Location)> {
        #[allow(unused_mut)] // Because we need mut in 2D but not in 3D.
        if let Some((mut proj, loc)) =
            CompositeShapeRef(self).project_local_point_and_get_location(point, max_dist, solid)
        {
            let seg_id = proj.subshape;

            // A point behind the outward pseudo-normal is inside.
            #[cfg(feature = "dim2")]
            if let Some(constraints) = self.segment_normal_constraints(seg_id) {
                let pseudo_normal = match loc {
                    SegmentPointLocation::OnVertex(i) => constraints.edges[i as usize],
                    SegmentPointLocation::OnEdge(_) => constraints.face,
                };
                proj.is_inside = (point - proj.point).dot(pseudo_normal) <= 0.0;

                if proj.is_inside && solid {
                    proj.point = point;
                }
            }

            Some((proj, (seg_id, loc)))
        } else {
            None
        }
    }
}

impl PointQueryWithLocation for TriMesh {
    type Location = (u32, TrianglePointLocation);

    #[inline]
    #[allow(unused_mut)] // Because we need mut in 3D but not in 2D.
    fn project_local_point_and_get_location(
        &self,
        point: Vector,
        solid: bool,
    ) -> (PointProjection, Self::Location) {
        self.project_local_point_and_get_location_with_max_dist(point, solid, Real::MAX)
            // No candidate: `point` (or `self`) isn’t finite. See
            // `Polyline::project_local_point_and_get_feature`.
            .unwrap_or((
                PointProjection::new(false, point),
                (0, TrianglePointLocation::OnVertex(0)),
            ))
    }

    /// Projects a point on `self`, with a maximum projection distance.
    fn project_local_point_and_get_location_with_max_dist(
        &self,
        point: Vector,
        solid: bool,
        max_dist: Real,
    ) -> Option<(PointProjection, Self::Location)> {
        #[allow(unused_mut)] // mut is needed in 3D.
        if let Some((mut proj, location)) =
            CompositeShapeRef(self).project_local_point_and_get_location(point, max_dist, solid)
        {
            let part_id = proj.subshape;

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
                    proj.is_inside = dpt.dot(pseudo_normal) <= 0.0;
                }
            }

            Some((proj, (part_id, location)))
        } else {
            None
        }
    }
}
