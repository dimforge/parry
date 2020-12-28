use crate::bounding_volume::SimdAABB;
use crate::math::{Point, Real, SimdReal, SIMD_WIDTH};
use crate::partitioning::{SimdBestFirstVisitStatus, SimdBestFirstVisitor};
use crate::query::{
    visitors::CompositePointContainmentTest, PointProjection, PointQuery, PointQueryWithLocation,
};
use crate::shape::{FeatureId, TriMesh, TrianglePointLocation};
use na;
use simba::simd::{SimdBool as _, SimdPartialOrd, SimdValue};

impl PointQuery for TriMesh {
    #[inline]
    fn project_local_point(&self, point: &Point<Real>, solid: bool) -> PointProjection {
        self.project_local_point_with_location(point, solid).0
    }

    #[inline]
    fn project_local_point_and_get_feature(
        &self,
        point: &Point<f32>,
    ) -> (PointProjection, FeatureId) {
        let (projection, (triangle_id, _location)) =
            self.project_local_point_with_location(point, false);

        /*
        let feature_id = match location {
            TrianglePointLocation::OnVertex(vid) => {
                FeatureId::Vertex(triangle_id)
            }
            TrianglePointLocation::OnEdge(eid, _) => {
                FeatureId::Edge(face.edges[triangle_local_id])
            }
            TrianglePointLocation::OnFace(_, _) => ,
            TrianglePointLocation::OnSolid => FeatureId::Unknown,
        };
         */
        // TODO: find a way to give proper feature ids to edges
        // and vertices.
        let feature_id = FeatureId::Face(triangle_id);
        (projection, feature_id)
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

impl PointQueryWithLocation for TriMesh {
    type Location = (u32, TrianglePointLocation);

    #[inline]
    fn project_local_point_with_location(
        &self,
        point: &Point<Real>,
        _: bool,
    ) -> (PointProjection, Self::Location) {
        let mut visitor = TriMeshPointProjVisitor {
            mesh: self,
            point,
            simd_point: Point::splat(*point),
        };

        self.quadtree().traverse_best_first(&mut visitor).unwrap().1
    }
}

/*
 * Visitors
 */
struct TriMeshPointProjVisitor<'a> {
    mesh: &'a TriMesh,
    point: &'a Point<Real>,
    simd_point: Point<SimdReal>,
}

impl<'a> SimdBestFirstVisitor<u32, SimdAABB> for TriMeshPointProjVisitor<'a> {
    type Result = (PointProjection, (u32, TrianglePointLocation));

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
                if (bitmask & (1 << ii)) != 0 && data[ii].is_some() {
                    let triangle = self.mesh.triangle(*data[ii].unwrap());
                    let (proj, extra_info) =
                        triangle.project_local_point_with_location(self.point, true);
                    results[ii] = Some((proj, (*data[ii].unwrap(), extra_info)));
                    weights[ii] = na::distance(self.point, &proj.point);
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
