use crate::bounding_volume::BoundingVolume;
use crate::math::{Pose, Real, Vector};
use crate::query::details::ShapeCastOptions;
use crate::query::{QueryDispatcher, Ray, ShapeCastHit, Unsupported};
use crate::shape::{HeightField, Shape};
#[cfg(feature = "dim3")]
use crate::{
    bounding_volume::Aabb,
    query::{PointQuery, RayCast},
    shape::{pseudo_normals::project_into_cone, HeightFieldCellStatus, HeightFieldFlags, Triangle},
};

/// Time Of Impact between a moving shape and a heightfield.
#[cfg(feature = "dim2")]
pub fn cast_shapes_heightfield_shape<D: ?Sized + QueryDispatcher>(
    dispatcher: &D,
    pos12: &Pose,
    vel12: Vector,
    heightfield1: &HeightField,
    g2: &dyn Shape,
    options: ShapeCastOptions,
) -> Result<Option<ShapeCastHit>, Unsupported> {
    let aabb2_1 = g2.compute_aabb(pos12).loosened(options.target_distance);
    let ray = Ray::new(aabb2_1.center(), vel12);

    let mut curr_range = heightfield1.unclamped_elements_range_in_local_aabb(&aabb2_1);
    // Enlarge the range by 1 to account for movement within a cell.
    let right = ray.dir.x > 0.0;

    if right {
        curr_range.end += 1;
    } else {
        curr_range.start -= 1;
    }

    let mut best_hit = None::<ShapeCastHit>;

    /*
     * Test the segment under the ray.
     */
    let clamped_curr_range = curr_range.start.clamp(0, heightfield1.num_cells() as isize) as usize
        ..curr_range.end.clamp(0, heightfield1.num_cells() as isize) as usize;
    for curr in clamped_curr_range {
        if let Some(seg) = heightfield1.segment_at(curr) {
            // TODO: pre-check using a ray-cast on the Aabbs first?
            if let Some(hit) = dispatcher.cast_shapes(pos12, vel12, &seg, g2, options)? {
                if hit.time_of_impact < best_hit.map(|h| h.time_of_impact).unwrap_or(Real::MAX) {
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
            curr_param = (cell_width * (curr_elt as Real) + start_x - ray.origin.x) / ray.dir.x;
        } else {
            curr_param = (ray.origin.x - cell_width * (curr_elt as Real) - start_x) / ray.dir.x;
            curr_elt -= 1;
        }

        if curr_param >= options.max_time_of_impact {
            break;
        }

        if let Some(seg) = heightfield1.segment_at(curr_elt as usize) {
            // TODO: pre-check using a ray-cast on the Aabbs first?
            if let Some(hit) = dispatcher.cast_shapes(pos12, vel12, &seg, g2, options)? {
                if hit.time_of_impact < best_hit.map(|h| h.time_of_impact).unwrap_or(Real::MAX) {
                    best_hit = Some(hit);
                }
            }
        }
    }

    Ok(best_hit)
}

/// Time Of Impact between a moving shape and a heightfield.
#[cfg(feature = "dim3")]
pub fn cast_shapes_heightfield_shape<D: ?Sized + QueryDispatcher>(
    dispatcher: &D,
    pos12: &Pose,
    vel12: Vector,
    heightfield1: &HeightField,
    g2: &dyn Shape,
    options: ShapeCastOptions,
) -> Result<Option<ShapeCastHit>, Unsupported> {
    let aabb1 = heightfield1.local_aabb();
    let mut aabb2_1 = g2.compute_aabb(pos12).loosened(options.target_distance);
    let ray = Ray::new(aabb2_1.center(), vel12);

    // Find the first hit between the aabbs.
    let hext2_1 = aabb2_1.half_extents();
    let msum = Aabb::new(aabb1.mins - hext2_1, aabb1.maxs + hext2_1);
    if let Some(time_of_impact) = msum.cast_local_ray(&ray, options.max_time_of_impact, true) {
        // Advance the aabb2 to the hit point.
        aabb2_1.mins += ray.dir * time_of_impact;
        aabb2_1.maxs += ray.dir * time_of_impact;
    } else {
        return Ok(None);
    }

    let (mut curr_range_i, mut curr_range_j) =
        heightfield1.unclamped_elements_range_in_local_aabb(&aabb2_1);
    let (ncells_i, ncells_j) = heightfield1.num_cells_ij();
    let mut best_hit = None::<ShapeCastHit>;

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

    let fix_internal_edges = heightfield1
        .flags()
        .contains(HeightFieldFlags::FIX_INTERNAL_EDGES);
    let mut hit_triangles = |i, j| {
        if i >= 0 && j >= 0 {
            let (tri_a, tri_b) = heightfield1.triangles_at(i as usize, j as usize);
            for tri in [tri_a, tri_b].into_iter().flatten() {
                // TODO: pre-check using a ray-cast on the Aabbs first?
                if let Some(mut hit) = dispatcher.cast_shapes(pos12, vel12, &tri, g2, options)? {
                    if hit.time_of_impact < best_hit.map(|h| h.time_of_impact).unwrap_or(Real::MAX)
                    {
                        // Each triangle is cast on its own, so the hit may lie on an internal edge,
                        // with a normal the heightfield's surface doesn't have.
                        if fix_internal_edges {
                            if let Some(normal1) = project_into_surface_cone(
                                heightfield1,
                                &tri,
                                hit.witness1,
                                hit.normal1,
                            ) {
                                hit.normal1 = normal1;
                                hit.normal2 = pos12.rotation.inverse() * -normal1;
                            }
                        }
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

    let mut cell = heightfield1.unclamped_cell_at_point(aabb2_1.center());

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

        if toi_x > options.max_time_of_impact && toi_z > options.max_time_of_impact {
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

/// `normal`, the normal of a hit at `witness` on `triangle` of `heightfield`, projected into the
/// normal cone of the surface there, as the triangles of `heightfield` within a band around
/// `witness` approximate it: a cone around the face normal reaching the normal of the triangle in
/// the band that lies farthest from it. Inside a face, on an edge to a coplanar triangle, or at a
/// vertex of a flat region, that leaves only the face normal; on a real crease, the normals
/// between its faces.
///
/// `None`, leaving `normal` as it is, where the band reaches the heightfield's border or a removed
/// triangle, as the cone opens towards the missing surface there, or if `normal` points to the
/// triangle's back.
#[cfg(feature = "dim3")]
fn project_into_surface_cone(
    heightfield: &HeightField,
    triangle: &Triangle,
    witness: Vector,
    mut normal: Vector,
) -> Option<Vector> {
    let face = triangle.normal()?;
    // The witness lies on the surface up to the cast's rounding. A wider band only widens the cone.
    let tolerance = 1.0e-3
        * triangle
            .edges()
            .iter()
            .map(|edge| edge.length())
            .fold(0.0, Real::max);
    let band = Aabb::new(
        witness - Vector::splat(tolerance),
        witness + Vector::splat(tolerance),
    );
    let (cells_i, cells_j) = heightfield.unclamped_elements_range_in_local_aabb(&band);

    let mut farthest = face;
    for i in cells_i {
        for j in cells_j.clone() {
            let (Ok(i), Ok(j)) = (usize::try_from(i), usize::try_from(j)) else {
                return None;
            };
            if i >= heightfield.nrows() || j >= heightfield.ncols() {
                return None;
            }
            let status = heightfield.cell_status(i, j);
            let (left, right) = heightfield.cell_triangles(i, j);
            for (other, removed) in [
                (left, HeightFieldCellStatus::LEFT_TRIANGLE_REMOVED),
                (right, HeightFieldCellStatus::RIGHT_TRIANGLE_REMOVED),
            ] {
                if status.contains(removed) {
                    if reaches_cylinder(&other, witness, tolerance) {
                        return None;
                    }
                    continue;
                }
                if other.distance_to_local_point(witness, true) <= tolerance {
                    let other_normal = other.normal()?;
                    if other_normal.dot(face) < farthest.dot(face) {
                        farthest = other_normal;
                    }
                }
            }
        }
    }

    let pseudo_normal = (face + farthest).try_normalize()?;
    project_into_cone(face, pseudo_normal, &mut normal).then_some(normal)
}

/// Whether `triangle` reaches the cylinder around `point` of radius and half-height `tolerance`,
/// which holds the ball of radius `tolerance` around it.
///
/// Only the triangle's part between the cylinder's heights is measured, seen from above, so heights
/// far from `point`, however extreme, cannot round it away. A non-finite height counts as reaching.
#[cfg(feature = "dim3")]
fn reaches_cylinder(triangle: &Triangle, point: Vector, tolerance: Real) -> bool {
    let corners = [triangle.a, triangle.b, triangle.c];
    if corners.iter().any(|corner| !corner.y.is_finite()) {
        return true;
    }
    let (low, high) = (point.y - tolerance, point.y + tolerance);

    // The triangle clipped to the slab between `low` and `high`: each edge adds its start, where
    // within the slab, and where it crosses the slab's bounds, in order along the edge.
    let mut part = [Vector::ZERO; 9];
    let mut len = 0;
    for k in 0..3 {
        let (a, b) = (corners[k], corners[(k + 1) % 3]);
        if (low..=high).contains(&a.y) {
            part[len] = a;
            len += 1;
        }
        let bounds = if a.y <= b.y { [low, high] } else { [high, low] };
        for bound in bounds {
            if (a.y < bound) != (b.y < bound) && a.y != bound && b.y != bound {
                // Halved, the heights' differences cannot overflow; only `x` and `z` are read below.
                let t = (bound * 0.5 - a.y * 0.5) / (b.y * 0.5 - a.y * 0.5);
                part[len] = Vector::new(a.x + (b.x - a.x) * t, bound, a.z + (b.z - a.z) * t);
                len += 1;
            }
        }
    }

    // The clipped part seen from above, relative to `point`: within `tolerance` of its boundary,
    // or with `point` inside it.
    let plan = |v: Vector| Vector::new(v.x - point.x, 0.0, v.z - point.z);
    let mut nearest = Real::MAX;
    let mut side = 0.0;
    let mut inside = true;
    for k in 0..len {
        let a = plan(part[k]);
        let edge = plan(part[(k + 1) % len]) - a;
        let length_squared = edge.length_squared();
        let t = if length_squared > 0.0 {
            (-a.dot(edge) / length_squared).clamp(0.0, 1.0)
        } else {
            0.0
        };
        nearest = nearest.min((a + edge * t).length());
        let cross = edge.z * a.x - edge.x * a.z;
        if cross != 0.0 {
            if side == 0.0 {
                side = cross.signum();
            } else if cross.signum() != side {
                inside = false;
            }
        }
    }
    len > 0 && (nearest <= tolerance || (inside && side != 0.0))
}

#[cfg(all(test, feature = "dim3"))]
mod test {
    use super::reaches_cylinder;
    use crate::math::{Real, Vector};
    use crate::shape::Triangle;

    /// The triangle of the cell between `x = 2` and `x = 3`, `z = -1` and `z = 0` where
    /// `x + z >= 2`, on the plane of `height`.
    fn removed(height: impl Fn(Real, Real) -> Real) -> Triangle {
        let corner = |x: Real, z: Real| Vector::new(x, height(x, z), z);
        Triangle::new(corner(3.0, -1.0), corner(3.0, 0.0), corner(2.0, 0.0))
    }

    #[test]
    fn a_triangle_reaches_the_cylinder_on_its_border() {
        let triangle = removed(|x, z| x + z);
        let on_border = Vector::new(2.5, 2.0, -0.5);
        assert!(reaches_cylinder(&triangle, on_border, 1.0e-3));
        let beside = Vector::new(2.45, 1.9, -0.55);
        assert!(!reaches_cylinder(&triangle, beside, 1.0e-3));
        let inside = Vector::new(2.8, 2.6, -0.2);
        assert!(reaches_cylinder(&triangle, inside, 1.0e-3));
    }

    #[test]
    fn a_far_corner_does_not_round_the_border_away() {
        let triangle = removed(|x, z| if x == 3.0 && z == 0.0 { 1.0e20 } else { x + z });
        let on_border = Vector::new(2.5, 2.0, -0.5);
        assert!(reaches_cylinder(&triangle, on_border, 1.0e-3));
        let below = removed(|x, z| if x == 3.0 && z == 0.0 { -1.0e20 } else { x + z });
        assert!(reaches_cylinder(&below, on_border, 1.0e-3));
        for height in [Real::INFINITY, Real::NAN] {
            let triangle = removed(|x, z| if x == 3.0 && z == 0.0 { height } else { x + z });
            assert!(reaches_cylinder(&triangle, on_border, 1.0e-3));
        }
    }

    #[test]
    fn heights_of_opposite_extremes_do_not_overflow_the_clipping() {
        // The point lies on the edge from 3e38 below to 3e38 above, whose rise overflows.
        let triangle = Triangle::new(
            Vector::new(3.0, -3.0e38, -1.0),
            Vector::new(3.0, 3.0e38, 0.0),
            Vector::new(2.0, 3.0e38, 0.0),
        );
        assert!(reaches_cylinder(
            &triangle,
            Vector::new(3.0, 0.0, -0.5),
            1.0e-3
        ));
    }

    #[test]
    fn a_steep_triangle_metres_away_does_not_reach_the_cylinder() {
        // Seen from above, the triangle lies 3.5 cm from the point; in space, 5 m.
        let triangle = removed(|x, z| 100.0 * (x + z));
        let point = Vector::new(2.0, 195.0, -0.05);
        assert!(!reaches_cylinder(&triangle, point, 0.1));
    }
}

/// Time Of Impact between a moving shape and a heightfield.
pub fn cast_shapes_shape_heightfield<D: ?Sized + QueryDispatcher>(
    dispatcher: &D,
    pos12: &Pose,
    vel12: Vector,
    g1: &dyn Shape,
    heightfield2: &HeightField,
    options: ShapeCastOptions,
) -> Result<Option<ShapeCastHit>, Unsupported> {
    Ok(cast_shapes_heightfield_shape(
        dispatcher,
        &pos12.inverse(),
        -(pos12.rotation.inverse() * vel12),
        heightfield2,
        g1,
        options,
    )?
    .map(|hit| hit.swapped()))
}
