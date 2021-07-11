use crate::bounding_volume::SimdAABB;
use crate::math::{Real, SimdBool, SimdReal, SIMD_WIDTH};
use crate::partitioning::{SimdBestFirstVisitStatus, SimdBestFirstVisitor};
use crate::query::{Ray, RayCast, RayIntersection, SimdRay};
use crate::shape::{Compound, FeatureId, Polyline, TriMesh, TypedSimdCompositeShape};
use simba::simd::{SimdBool as _, SimdPartialOrd, SimdValue};

impl RayCast for TriMesh {
    #[inline]
    fn cast_local_ray(&self, ray: &Ray, max_toi: Real, solid: bool) -> Option<Real> {
        let mut visitor = RayCompositeShapeToiBestFirstVisitor::new(self, ray, max_toi, solid);

        self.qbvh()
            .traverse_best_first(&mut visitor)
            .map(|res| res.1 .1)
    }

    #[inline]
    fn cast_local_ray_and_get_normal(
        &self,
        ray: &Ray,
        max_toi: Real,
        solid: bool,
    ) -> Option<RayIntersection> {
        let mut visitor =
            RayCompositeShapeToiAndNormalBestFirstVisitor::new(self, ray, max_toi, solid);

        self.qbvh()
            .traverse_best_first(&mut visitor)
            .map(|(_, (best, mut res))| {
                // We hit a backface.
                // NOTE: we need this for `TriMesh::is_backface` to work properly.
                if res.feature == FeatureId::Face(1) {
                    res.feature = FeatureId::Face(best + self.indices().len() as u32)
                } else {
                    res.feature = FeatureId::Face(best);
                }
                res
            })
    }
}

impl RayCast for Polyline {
    #[inline]
    fn cast_local_ray(&self, ray: &Ray, max_toi: Real, solid: bool) -> Option<Real> {
        let mut visitor = RayCompositeShapeToiBestFirstVisitor::new(self, ray, max_toi, solid);

        self.qbvh()
            .traverse_best_first(&mut visitor)
            .map(|res| res.1 .1)
    }

    #[inline]
    fn cast_local_ray_and_get_normal(
        &self,
        ray: &Ray,
        max_toi: Real,
        solid: bool,
    ) -> Option<RayIntersection> {
        let mut visitor =
            RayCompositeShapeToiAndNormalBestFirstVisitor::new(self, ray, max_toi, solid);

        self.qbvh()
            .traverse_best_first(&mut visitor)
            .map(|(_, (_, res))| res)
    }
}

impl RayCast for Compound {
    #[inline]
    fn cast_local_ray(&self, ray: &Ray, max_toi: Real, solid: bool) -> Option<Real> {
        let mut visitor = RayCompositeShapeToiBestFirstVisitor::new(self, ray, max_toi, solid);

        self.qbvh()
            .traverse_best_first(&mut visitor)
            .map(|res| res.1 .1)
    }

    #[inline]
    fn cast_local_ray_and_get_normal(
        &self,
        ray: &Ray,
        max_toi: Real,
        solid: bool,
    ) -> Option<RayIntersection> {
        let mut visitor =
            RayCompositeShapeToiAndNormalBestFirstVisitor::new(self, ray, max_toi, solid);

        self.qbvh()
            .traverse_best_first(&mut visitor)
            .map(|(_, (_, res))| res)
    }
}

/*
 * Visitors
 */
/// A visitor for casting a ray on a composite shape.
pub struct RayCompositeShapeToiBestFirstVisitor<'a, S> {
    shape: &'a S,
    ray: &'a Ray,
    simd_ray: SimdRay,
    max_toi: Real,
    solid: bool,
}

impl<'a, S> RayCompositeShapeToiBestFirstVisitor<'a, S> {
    /// Initialize a visitor for casting a ray on a composite shape.
    pub fn new(shape: &'a S, ray: &'a Ray, max_toi: Real, solid: bool) -> Self {
        Self {
            shape,
            ray,
            simd_ray: SimdRay::splat(*ray),
            max_toi,
            solid,
        }
    }
}

impl<'a, S> SimdBestFirstVisitor<S::PartId, SimdAABB>
    for RayCompositeShapeToiBestFirstVisitor<'a, S>
where
    S: TypedSimdCompositeShape,
{
    type Result = (S::PartId, Real);

    #[inline]
    fn visit(
        &mut self,
        best: Real,
        aabb: &SimdAABB,
        data: Option<[Option<&S::PartId>; SIMD_WIDTH]>,
    ) -> SimdBestFirstVisitStatus<Self::Result> {
        let (hit, toi) = aabb.cast_local_ray(&self.simd_ray, SimdReal::splat(self.max_toi));

        if let Some(data) = data {
            let mut weights = [0.0; SIMD_WIDTH];
            let mut mask = [false; SIMD_WIDTH];
            let mut results = [None; SIMD_WIDTH];

            let better_toi = toi.simd_lt(SimdReal::splat(best));
            let bitmask = (hit & better_toi).bitmask();

            for ii in 0..SIMD_WIDTH {
                if (bitmask & (1 << ii)) != 0 && data[ii].is_some() {
                    let part_id = *data[ii].unwrap();
                    self.shape
                        .map_typed_part_at(part_id, |part_pos, part_shape| {
                            let toi = if let Some(part_pos) = part_pos {
                                part_shape.cast_ray(part_pos, &self.ray, self.max_toi, self.solid)
                            } else {
                                part_shape.cast_local_ray(&self.ray, self.max_toi, self.solid)
                            };
                            if let Some(toi) = toi {
                                results[ii] = Some((part_id, toi));
                                mask[ii] = true;
                                weights[ii] = toi;
                            }
                        })
                }
            }

            SimdBestFirstVisitStatus::MaybeContinue {
                weights: SimdReal::from(weights),
                mask: SimdBool::from(mask),
                results,
            }
        } else {
            SimdBestFirstVisitStatus::MaybeContinue {
                weights: toi,
                mask: hit,
                results: [None; SIMD_WIDTH],
            }
        }
    }
}

/// A visitor for casting a ray on a composite shape.
pub struct RayCompositeShapeToiAndNormalBestFirstVisitor<'a, S> {
    shape: &'a S,
    ray: &'a Ray,
    simd_ray: SimdRay,
    max_toi: Real,
    solid: bool,
}

impl<'a, S> RayCompositeShapeToiAndNormalBestFirstVisitor<'a, S> {
    /// Initialize a visitor for casting a ray on a composite shape.
    pub fn new(shape: &'a S, ray: &'a Ray, max_toi: Real, solid: bool) -> Self {
        Self {
            shape,
            ray,
            simd_ray: SimdRay::splat(*ray),
            max_toi,
            solid,
        }
    }
}

impl<'a, S> SimdBestFirstVisitor<S::PartId, SimdAABB>
    for RayCompositeShapeToiAndNormalBestFirstVisitor<'a, S>
where
    S: TypedSimdCompositeShape,
{
    type Result = (S::PartId, RayIntersection);

    #[inline]
    fn visit(
        &mut self,
        best: Real,
        aabb: &SimdAABB,
        data: Option<[Option<&S::PartId>; SIMD_WIDTH]>,
    ) -> SimdBestFirstVisitStatus<Self::Result> {
        let (hit, toi) = aabb.cast_local_ray(&self.simd_ray, SimdReal::splat(self.max_toi));

        if let Some(data) = data {
            let mut weights = [0.0; SIMD_WIDTH];
            let mut mask = [false; SIMD_WIDTH];
            let mut results = [None; SIMD_WIDTH];

            let better_toi = toi.simd_lt(SimdReal::splat(best));
            let bitmask = (hit & better_toi).bitmask();

            for ii in 0..SIMD_WIDTH {
                if (bitmask & (1 << ii)) != 0 && data[ii].is_some() {
                    self.shape
                        .map_typed_part_at(*data[ii].unwrap(), |part_pos, part_shape| {
                            let result = if let Some(part_pos) = part_pos {
                                part_shape.cast_ray_and_get_normal(
                                    part_pos,
                                    &self.ray,
                                    self.max_toi,
                                    self.solid,
                                )
                            } else {
                                part_shape.cast_local_ray_and_get_normal(
                                    &self.ray,
                                    self.max_toi,
                                    self.solid,
                                )
                            };

                            if let Some(result) = result {
                                results[ii] = Some((*data[ii].unwrap(), result));
                                mask[ii] = true;
                                weights[ii] = result.toi;
                            }
                        });
                }
            }

            SimdBestFirstVisitStatus::MaybeContinue {
                weights: SimdReal::from(weights),
                mask: SimdBool::from(mask),
                results,
            }
        } else {
            SimdBestFirstVisitStatus::MaybeContinue {
                weights: toi,
                mask: hit,
                results: [None; SIMD_WIDTH],
            }
        }
    }
}
