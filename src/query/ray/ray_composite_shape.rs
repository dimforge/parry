use crate::bounding_volume::SimdAABB;
use crate::math::{Real, SimdBool, SimdReal, SIMD_WIDTH};
use crate::partitioning::{SimdBestFirstVisitStatus, SimdBestFirstVisitor};
use crate::query::{Ray, RayCast, RayIntersection, SimdRay};
use crate::shape::{Compound, FeatureId, Polyline, TriMesh, TypedSimdCompositeShape};
use simba::simd::{SimdBool as _, SimdPartialOrd, SimdValue};

impl RayCast for TriMesh {
    #[inline]
    fn cast_local_ray(&self, ray: &Ray, max_toi: Real, solid: bool) -> Option<Real> {
        let mut visitor = CompositeShapeRayToiVisitor {
            shape: self,
            ray,
            simd_ray: SimdRay::splat(*ray),
            max_toi,
            solid,
        };

        self.quadtree()
            .traverse_best_first(&mut visitor)
            .map(|res| res.1)
    }

    #[inline]
    fn cast_local_ray_and_get_normal(
        &self,
        ray: &Ray,
        max_toi: Real,
        solid: bool,
    ) -> Option<RayIntersection> {
        let mut visitor = CompositeShapeRayToiAndNormalVisitor {
            shape: self,
            ray,
            simd_ray: SimdRay::splat(*ray),
            max_toi,
            solid,
        };

        self.quadtree()
            .traverse_best_first(&mut visitor)
            .map(|(_, (best, mut res))| {
                res.feature = FeatureId::Face(best);
                res
            })
    }
}

impl RayCast for Polyline {
    #[inline]
    fn cast_local_ray(&self, ray: &Ray, max_toi: Real, solid: bool) -> Option<Real> {
        let mut visitor = CompositeShapeRayToiVisitor {
            shape: self,
            ray,
            simd_ray: SimdRay::splat(*ray),
            max_toi,
            solid,
        };

        self.quadtree()
            .traverse_best_first(&mut visitor)
            .map(|res| res.1)
    }

    #[inline]
    fn cast_local_ray_and_get_normal(
        &self,
        ray: &Ray,
        max_toi: Real,
        solid: bool,
    ) -> Option<RayIntersection> {
        let mut visitor = CompositeShapeRayToiAndNormalVisitor {
            shape: self,
            ray,
            simd_ray: SimdRay::splat(*ray),
            max_toi,
            solid,
        };

        self.quadtree()
            .traverse_best_first(&mut visitor)
            .map(|(_, (_, res))| res)
    }
}

impl RayCast for Compound {
    #[inline]
    fn cast_local_ray(&self, ray: &Ray, max_toi: Real, solid: bool) -> Option<Real> {
        let mut visitor = CompositeShapeRayToiVisitor {
            shape: self,
            ray,
            simd_ray: SimdRay::splat(*ray),
            max_toi,
            solid,
        };

        self.quadtree()
            .traverse_best_first(&mut visitor)
            .map(|res| res.1)
    }

    #[inline]
    fn cast_local_ray_and_get_normal(
        &self,
        ray: &Ray,
        max_toi: Real,
        solid: bool,
    ) -> Option<RayIntersection> {
        let mut visitor = CompositeShapeRayToiAndNormalVisitor {
            shape: self,
            ray,
            simd_ray: SimdRay::splat(*ray),
            max_toi,
            solid,
        };

        self.quadtree()
            .traverse_best_first(&mut visitor)
            .map(|(_, (_, res))| res)
    }
}

/*
 * Visitors
 */
struct CompositeShapeRayToiVisitor<'a, S> {
    shape: &'a S,
    ray: &'a Ray,
    simd_ray: SimdRay,
    max_toi: Real,
    solid: bool,
}

impl<'a, S> SimdBestFirstVisitor<u32, SimdAABB> for CompositeShapeRayToiVisitor<'a, S>
where
    S: TypedSimdCompositeShape,
{
    type Result = Real;

    #[inline]
    fn visit(
        &mut self,
        best: Real,
        aabb: &SimdAABB,
        data: Option<[Option<&u32>; SIMD_WIDTH]>,
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
                            let toi = if let Some(part_pos) = part_pos {
                                part_shape.cast_ray(part_pos, &self.ray, self.max_toi, self.solid)
                            } else {
                                part_shape.cast_local_ray(&self.ray, self.max_toi, self.solid)
                            };
                            if let Some(toi) = toi {
                                results[ii] = Some(toi);
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

struct CompositeShapeRayToiAndNormalVisitor<'a, S> {
    shape: &'a S,
    ray: &'a Ray,
    simd_ray: SimdRay,
    max_toi: Real,
    solid: bool,
}

impl<'a, S> SimdBestFirstVisitor<u32, SimdAABB> for CompositeShapeRayToiAndNormalVisitor<'a, S>
where
    S: TypedSimdCompositeShape,
{
    type Result = (u32, RayIntersection);

    #[inline]
    fn visit(
        &mut self,
        best: Real,
        aabb: &SimdAABB,
        data: Option<[Option<&u32>; SIMD_WIDTH]>,
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
