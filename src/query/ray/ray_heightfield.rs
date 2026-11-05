use crate::math::Real;
#[cfg(feature = "dim2")]
use crate::query;
use crate::query::{Ray, RayCast, RayIntersection};
#[cfg(feature = "dim2")]
use crate::shape::FeatureId;
use crate::shape::HeightField;

#[cfg(feature = "dim2")]
impl RayCast for HeightField {
    #[inline]
    fn cast_local_ray_and_get_normal(
        &self,
        ray: &Ray,
        max_time_of_impact: Real,
        _: bool,
    ) -> Option<RayIntersection> {
        let aabb = self.local_aabb();
        let (min_t, mut max_t) = aabb.clip_ray_parameters(ray)?;

        if min_t > max_time_of_impact {
            return None;
        }

        max_t = max_t.min(max_time_of_impact);

        let clip_ray_a = ray.point_at(min_t);

        // None may happen due to slight numerical errors.
        let mut curr = self.cell_at_point(clip_ray_a).unwrap_or_else(|| {
            if ray.origin.x > 0.0 {
                self.num_cells() - 1
            } else {
                0_usize
            }
        });

        /*
         * Test the segment under the ray.
         */
        let cast_first_segment = |i: usize| {
            let seg = self.segment_at(i)?;
            let (s, t) = if ray.dir.x == 0.0 {
                // With a positive scale, `cell_at_point` puts a vertical ray within the segment's
                // span, so this `t` lies in `[0, 1]`. The line-line parameters can round it out of
                // that range near the segment's ends, on a sloped segment.
                let t = (ray.origin.x - seg.a.x) / (seg.b.x - seg.a.x);
                let y = seg.a.y + (seg.b.y - seg.a.y) * t;
                ((y - ray.origin.y) / ray.dir.y, t)
            } else {
                query::details::closest_points_line_line_parameters(
                    ray.origin,
                    ray.dir,
                    seg.a,
                    seg.scaled_direction(),
                )
            };
            if s >= 0.0 && t >= 0.0 && t <= 1.0 {
                let n = seg.normal().unwrap();
                // The feature is the hit segment's own side; `subshape` says which segment.
                let fid = if n.dot(ray.dir) > 0.0 { 1 } else { 0 };
                Some(RayIntersection::new(s, n, FeatureId::Face(fid)).with_subshape(i as u32))
            } else {
                None
            }
        };

        if let Some(inter) = cast_first_segment(curr) {
            // Cast succeeded on the first element!
            return Some(inter);
        }

        /*
         * Test other segments in the path of the ray.
         */
        if ray.dir.x == 0.0 {
            // A vertical ray on a grid point touches the segments on both of its sides, and the
            // one `cell_at_point` gives may be removed.
            return self
                .cell_beside(curr, ray.origin.x)
                .and_then(cast_first_segment);
        }

        let right = ray.dir.x > 0.0;
        let cell_width = self.cell_width();
        let start_x = self.start_x();

        while (right && curr < self.num_cells()) || (!right && curr > 0) {
            let curr_param;

            if right {
                curr += 1;
                curr_param = (cell_width * (curr as Real) + start_x - ray.origin.x) / ray.dir.x;
            } else {
                curr_param = (ray.origin.x - cell_width * (curr as Real) - start_x) / ray.dir.x;
                curr -= 1;
            }

            if curr_param >= max_t {
                // The part of the ray after max_t is outside of the heightfield Aabb.
                return None;
            }

            if let Some(seg) = self.segment_at(curr) {
                // TODO: test the y-coordinates (equivalent to an Aabb test) before actually computing the intersection.
                let (s, t) = query::details::closest_points_line_line_parameters(
                    ray.origin,
                    ray.dir,
                    seg.a,
                    seg.scaled_direction(),
                );

                if t >= 0.0 && t <= 1.0 && s <= max_time_of_impact {
                    let n = seg.normal().unwrap();
                    // The feature is the hit segment's own side; `subshape` says which segment.
                    let fid = if n.dot(ray.dir) > 0.0 { 1 } else { 0 };
                    return Some(
                        RayIntersection::new(s, n, FeatureId::Face(fid)).with_subshape(curr as u32),
                    );
                }
            }
        }

        None
    }
}

#[cfg(feature = "dim3")]
impl RayCast for HeightField {
    #[inline]
    fn cast_local_ray_and_get_normal(
        &self,
        ray: &Ray,
        max_time_of_impact: Real,
        solid: bool,
    ) -> Option<RayIntersection> {
        use num_traits::Bounded;

        let aabb = self.local_aabb();
        let (min_t, mut max_t) = aabb.clip_ray_parameters(ray)?;
        max_t = max_t.min(max_time_of_impact);
        let clip_ray_a = ray.point_at(min_t);
        let mut cell = match self.cell_at_point(clip_ray_a) {
            Some(cell) => cell,
            // None may happen due to slight numerical errors.
            None => {
                let i = if ray.origin.z > 0.0 {
                    self.nrows() - 1
                } else {
                    0
                };

                let j = if ray.origin.x > 0.0 {
                    self.ncols() - 1
                } else {
                    0
                };

                (i, j)
            }
        };

        // The nearest hit among the triangles of the cell `(i, j)`.
        let cast_cell = |i: usize, j: usize| {
            let tris = self.triangles_at(i, j);
            let inter1 = tris
                .0
                .and_then(|tri| tri.cast_local_ray_and_get_normal(ray, max_time_of_impact, solid));
            let inter2 = tris
                .1
                .and_then(|tri| tri.cast_local_ray_and_get_normal(ray, max_time_of_impact, solid));

            // The feature stays the one the triangle reported; `subshape` says which triangle it
            // belongs to.
            match (inter1, inter2) {
                (Some(mut inter1), Some(mut inter2)) => {
                    if inter1.time_of_impact < inter2.time_of_impact {
                        inter1.subshape = self.triangle_id(i, j, true);
                        Some(inter1)
                    } else {
                        inter2.subshape = self.triangle_id(i, j, false);
                        Some(inter2)
                    }
                }
                (Some(mut inter), None) => {
                    inter.subshape = self.triangle_id(i, j, true);
                    Some(inter)
                }
                (None, Some(mut inter)) => {
                    inter.subshape = self.triangle_id(i, j, false);
                    Some(inter)
                }
                (None, None) => None,
            }
        };

        // A ray that runs exactly along a grid line, never crossing it, touches the cells on both
        // of its sides: the triangles' rounding decides which of the two reports a hit on their
        // shared edge, and the cell `cell_at_point` gives may be removed. Both meet the ray only
        // on that edge, so either one's hit will do. The ray stays on that line as it walks, so
        // the column (row) beside it stays the same.
        let column_beside = if ray.dir.x == 0.0 {
            self.column_beside(cell.1, clip_ray_a.x)
        } else {
            None
        };
        let row_beside = if ray.dir.z == 0.0 {
            self.row_beside(cell.0, clip_ray_a.z)
        } else {
            None
        };

        loop {
            for i in core::iter::once(cell.0).chain(row_beside) {
                for j in core::iter::once(cell.1).chain(column_beside) {
                    if let Some(inter) = cast_cell(i, j) {
                        return Some(inter);
                    }
                }
            }

            /*
             * Find the next cell to cast the ray on.
             */
            let (toi_x, right) = if ray.dir.x > 0.0 {
                let x = self.x_at(cell.1 + 1);
                ((x - ray.origin.x) / ray.dir.x, true)
            } else if ray.dir.x < 0.0 {
                let x = self.x_at(cell.1);
                ((x - ray.origin.x) / ray.dir.x, false)
            } else {
                (Real::max_value(), false)
            };

            let (toi_z, down) = if ray.dir.z > 0.0 {
                let z = self.z_at(cell.0 + 1);
                ((z - ray.origin.z) / ray.dir.z, true)
            } else if ray.dir.z < 0.0 {
                let z = self.z_at(cell.0);
                ((z - ray.origin.z) / ray.dir.z, false)
            } else {
                (Real::max_value(), false)
            };

            if toi_x > max_t && toi_z > max_t {
                break;
            }

            if toi_x >= 0.0 && toi_x < toi_z {
                if right {
                    cell.1 += 1
                } else if cell.1 > 0 {
                    cell.1 -= 1
                } else {
                    break;
                }
            } else if toi_z >= 0.0 {
                if down {
                    cell.0 += 1
                } else if cell.0 > 0 {
                    cell.0 -= 1
                } else {
                    break;
                }
            } else {
                break;
            }

            if cell.0 >= self.nrows() || cell.1 >= self.ncols() {
                break;
            }
        }

        None
    }
}
