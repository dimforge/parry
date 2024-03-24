use crate::math::{Isometry, Real, Vector};
use crate::query::{QueryDispatcher, Ray, Unsupported, TOI};
use crate::shape::{GenericHeightField, HeightFieldStorage, Shape};
#[cfg(feature = "dim3")]
use crate::{bounding_volume::Aabb, query::RayCast};

/// Time Of Impact between a moving shape and a heightfield.
#[cfg(feature = "dim2")]
pub fn time_of_impact_heightfield_shape<Storage, D: ?Sized>(
    dispatcher: &D,
    pos12: &Isometry<Real>,
    vel12: &Vector<Real>,
    heightfield1: &GenericHeightField<Storage>,
    g2: &dyn Shape,
    max_toi: Real,
    stop_at_penetration: bool,
) -> Result<Option<TOI>, Unsupported>
where
    Storage: HeightFieldStorage,
    D: QueryDispatcher,
{
    let aabb2_1 = g2.compute_aabb(pos12);
    let ray = Ray::new(aabb2_1.center(), *vel12);

    let mut curr_range = heightfield1.unclamped_elements_range_in_local_aabb(&aabb2_1);
    // Enlarge the range by 1 to account for movement within a cell.
    let right = ray.dir.x > 0.0;

    if right {
        curr_range.end += 1;
    } else {
        curr_range.start -= 1;
    }

    let mut best_hit = None::<TOI>;

    /*
     * Test the segment under the ray.
     */
    let clamped_curr_range = curr_range.start.clamp(0, heightfield1.num_cells() as isize) as usize
        ..curr_range.end.clamp(0, heightfield1.num_cells() as isize) as usize;
    for curr in clamped_curr_range {
        if let Some(seg) = heightfield1.segment_at(curr) {
            // TODO: pre-check using a ray-cast on the Aabbs first?
            if let Some(hit) =
                dispatcher.time_of_impact(pos12, vel12, &seg, g2, max_toi, stop_at_penetration)?
            {
                if hit.toi < best_hit.map(|toi| toi.toi).unwrap_or(Real::MAX) {
                    best_hit = Some(hit);
                }
            }
        }
    }

    /*
     * Test other segments in the path of the ray.
     */
    if ray.dir.x == 0.0 {
        return Ok(best_hit);
    }

    let cell_width = heightfield1.cell_width();
    let start_x = heightfield1.start_x();

    let mut curr_elt = if right {
        (curr_range.end - 1).max(0)
    } else {
        curr_range.start.min(heightfield1.num_cells() as isize - 1)
    };

    while (right && curr_elt < heightfield1.num_cells() as isize - 1) || (!right && curr_elt > 0) {
        let curr_param;

        if right {
            curr_elt += 1;
            curr_param = (cell_width * na::convert::<f64, Real>(curr_elt as f64) + start_x
                - ray.origin.x)
                / ray.dir.x;
        } else {
            curr_param =
                (ray.origin.x - cell_width * na::convert::<f64, Real>(curr_elt as f64) - start_x)
                    / ray.dir.x;
            curr_elt -= 1;
        }

        if curr_param >= max_toi {
            break;
        }

        if let Some(seg) = heightfield1.segment_at(curr_elt as usize) {
            // TODO: pre-check using a ray-cast on the Aabbs first?
            if let Some(hit) =
                dispatcher.time_of_impact(pos12, vel12, &seg, g2, max_toi, stop_at_penetration)?
            {
                if hit.toi < best_hit.map(|toi| toi.toi).unwrap_or(Real::MAX) {
                    best_hit = Some(hit);
                }
            }
        }
    }

    Ok(best_hit)
}

/// Time Of Impact between a moving shape and a heightfield.
#[cfg(feature = "dim3")]
pub fn time_of_impact_heightfield_shape<Storage, D: ?Sized>(
    dispatcher: &D,
    pos12: &Isometry<Real>,
    vel12: &Vector<Real>,
    heightfield1: &GenericHeightField<Storage>,
    g2: &dyn Shape,
    max_toi: Real,
    stop_at_penetration: bool,
) -> Result<Option<TOI>, Unsupported>
where
    Storage: HeightFieldStorage,
    D: QueryDispatcher,
{
    let aabb1 = heightfield1.local_aabb();
    let mut aabb2_1 = g2.compute_aabb(pos12);
    let ray = Ray::new(aabb2_1.center(), *vel12);

    // Find the first hit between the aabbs.
    let hext2_1 = aabb2_1.half_extents();
    let msum = Aabb::new(aabb1.mins - hext2_1, aabb1.maxs + hext2_1);
    if let Some(toi) = msum.cast_local_ray(&ray, max_toi, true) {
        // Advance the aabb2 to the hit point.
        aabb2_1.mins += ray.dir * toi;
        aabb2_1.maxs += ray.dir * toi;
    } else {
        return Ok(None);
    }

    let (mut curr_range_i, mut curr_range_j) =
        heightfield1.unclamped_elements_range_in_local_aabb(&aabb2_1);
    let (ncells_i, ncells_j) = heightfield1.num_cells_ij();
    let mut best_hit = None::<TOI>;

    /*
     * Enlarge the ranges by 1 to account for any movement within one cell.
     */
    if ray.dir.z > 0.0 {
        curr_range_i.end += 1;
    } else if ray.dir.z < 0.0 {
        curr_range_i.start -= 1;
    }

    if ray.dir.x > 0.0 {
        curr_range_j.end += 1;
    } else if ray.dir.x < 0.0 {
        curr_range_j.start -= 1;
    }

    /*
     * Test the segment under the ray.
     */
    let clamped_curr_range_i = curr_range_i.start.clamp(0, ncells_i as isize)
        ..curr_range_i.end.clamp(0, ncells_i as isize);
    let clamped_curr_range_j = curr_range_j.start.clamp(0, ncells_j as isize)
        ..curr_range_j.end.clamp(0, ncells_j as isize);

    let mut hit_triangles = |i, j| {
        if i >= 0 && j >= 0 {
            let (tri_a, tri_b) = heightfield1.triangles_at(i as usize, j as usize);
            for tri in [tri_a, tri_b].into_iter().flatten() {
                // TODO: pre-check using a ray-cast on the Aabbs first?
                if let Some(hit) = dispatcher.time_of_impact(
                    pos12,
                    vel12,
                    &tri,
                    g2,
                    max_toi,
                    stop_at_penetration,
                )? {
                    if hit.toi < best_hit.map(|toi| toi.toi).unwrap_or(Real::MAX) {
                        best_hit = Some(hit);
                    }
                }
            }
        }

        Ok(())
    };

    for i in clamped_curr_range_i {
        for j in clamped_curr_range_j.clone() {
            hit_triangles(i, j)?;
        }
    }

    if ray.dir.y == 0.0 {
        return Ok(best_hit);
    }

    let mut cell = heightfield1.unclamped_cell_at_point(&aabb2_1.center());

    loop {
        let prev_cell = cell;

        /*
         * Find the next cell to cast the ray on.
         */
        let toi_x = if ray.dir.x > 0.0 {
            let x = heightfield1.signed_x_at(cell.1 + 1);
            (x - ray.origin.x) / ray.dir.x
        } else if ray.dir.x < 0.0 {
            let x = heightfield1.signed_x_at(cell.1);
            (x - ray.origin.x) / ray.dir.x
        } else {
            Real::MAX
        };

        let toi_z = if ray.dir.z > 0.0 {
            let z = heightfield1.signed_z_at(cell.0 + 1);
            (z - ray.origin.z) / ray.dir.z
        } else if ray.dir.z < 0.0 {
            let z = heightfield1.signed_z_at(cell.0);
            (z - ray.origin.z) / ray.dir.z
        } else {
            Real::MAX
        };

        if toi_x > max_toi && toi_z > max_toi {
            break;
        }

        if toi_x >= 0.0 && toi_x <= toi_z {
            cell.1 += ray.dir.x.signum() as isize;
        }

        if toi_z >= 0.0 && toi_z <= toi_x {
            cell.0 += ray.dir.z.signum() as isize;
        }

        if cell == prev_cell {
            break;
        }

        let cell_diff = (cell.0 - prev_cell.0, cell.1 - prev_cell.1);
        curr_range_i.start += cell_diff.0;
        curr_range_i.end += cell_diff.0;
        curr_range_j.start += cell_diff.1;
        curr_range_j.end += cell_diff.1;

        let new_line_i = if cell_diff.0 > 0 {
            curr_range_i.end
        } else {
            curr_range_i.start
        };

        let new_line_j = if cell_diff.1 > 0 {
            curr_range_j.end
        } else {
            curr_range_j.start
        };

        let ignore_line_i = new_line_i < 0 || new_line_i >= ncells_i as isize;
        let ignore_line_j = new_line_j < 0 || new_line_j >= ncells_j as isize;

        if ignore_line_i && ignore_line_j {
            break;
        }

        if !ignore_line_i && cell_diff.0 != 0 {
            for j in curr_range_j.clone() {
                hit_triangles(new_line_i, j)?;
            }
        }

        if !ignore_line_j && cell_diff.1 != 0 {
            for i in curr_range_i.clone() {
                hit_triangles(i, new_line_j)?;
            }
        }
    }

    Ok(best_hit)
}

/// Time Of Impact between a moving shape and a heightfield.
pub fn time_of_impact_shape_heightfield<Storage, D: ?Sized>(
    dispatcher: &D,
    pos12: &Isometry<Real>,
    vel12: &Vector<Real>,
    g1: &dyn Shape,
    heightfield2: &GenericHeightField<Storage>,
    max_toi: Real,
    stop_at_penetration: bool,
) -> Result<Option<TOI>, Unsupported>
where
    Storage: HeightFieldStorage,
    D: QueryDispatcher,
{
    Ok(time_of_impact_heightfield_shape(
        dispatcher,
        &pos12.inverse(),
        &-pos12.inverse_transform_vector(vel12),
        heightfield2,
        g1,
        max_toi,
        stop_at_penetration,
    )?
    .map(|toi| toi.swapped()))
}
