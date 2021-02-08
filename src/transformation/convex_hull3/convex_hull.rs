use super::InitialMesh;
use super::TriangleFacet;
use crate::math::Real;
use crate::transformation::convex_hull_utils::indexed_support_point_nth;
use crate::transformation::convex_hull_utils::{indexed_support_point_id, normalize};
use crate::utils;
use na::{self, Point3};

/// Computes the convex hull of a set of 3d points.
pub fn convex_hull(points: &[Point3<Real>]) -> (Vec<Point3<Real>>, Vec<[u32; 3]>) {
    if points.is_empty() {
        return (Vec::new(), Vec::new());
    }

    // print_buildable_vec("input", points);

    let mut normalized_points = points.to_vec();
    let _ = normalize(&mut normalized_points[..]);

    let mut undecidable_points = Vec::new();
    let mut silhouette_loop_facets_and_idx = Vec::new();
    let mut removed_facets = Vec::new();

    let mut triangles;

    match super::get_initial_mesh(points, &mut normalized_points[..], &mut undecidable_points) {
        InitialMesh::Facets(facets) => {
            triangles = facets;
        }
        InitialMesh::ResultMesh(vertices, indices) => {
            return (vertices, indices);
        }
    }

    let mut i = 0;
    while i != triangles.len() {
        silhouette_loop_facets_and_idx.clear();

        if !triangles[i].valid || triangles[i].affinely_dependent {
            i = i + 1;
            continue;
        }

        // FIXME: use triangles[i].furthest_point instead.
        let pt_id = indexed_support_point_id(
            &triangles[i].normal,
            &normalized_points[..],
            triangles[i].visible_points[..].iter().copied(),
        );

        if let Some(point) = pt_id {
            triangles[i].valid = false;

            removed_facets.clear();
            removed_facets.push(i);

            for j in 0usize..3 {
                // println!(">> loop;");
                compute_silhouette(
                    triangles[i].adj[j],
                    triangles[i].indirect_adj_id[j],
                    point,
                    &mut silhouette_loop_facets_and_idx,
                    &normalized_points[..],
                    &mut removed_facets,
                    &mut triangles[..],
                );
            }

            // In some degenerate cases (because of float rounding problems), the silhouette may:
            // 1. Contain self-intersections (i.e. a single vertex is used by more than two edges).
            // 2. Contain multiple disjoint (but nested) loops.
            fix_silhouette_topology(
                &normalized_points,
                &mut silhouette_loop_facets_and_idx,
                &mut removed_facets,
                &mut triangles[..],
            );

            // Check that the silhouette is valid.
            // FIXME: remove this debug code.
            // {
            //     for (facet, id) in &silhouette_loop_facets_and_idx {
            //         assert!(triangles[*facet].valid);
            //         assert!(!triangles[triangles[*facet].adj[*id]].valid);
            //     }
            // }

            if silhouette_loop_facets_and_idx.is_empty() {
                // Due to inaccuracies, the silhouette could not be computed
                // (the point seems to be visible from… every triangle).
                let mut any_valid = false;
                for j in i + 1..triangles.len() {
                    if triangles[j].valid && !triangles[j].affinely_dependent {
                        any_valid = true;
                    }
                }

                if any_valid {
                    panic!(
                        "Warning: exitting an unfinished work: {}, {}",
                        normalized_points[point],
                        triangles.len()
                    );
                }

                // FIXME: this is very harsh.
                triangles[i].valid = true;
                break;
            }

            attach_and_push_facets(
                &silhouette_loop_facets_and_idx[..],
                point,
                &normalized_points[..],
                &mut triangles,
                &removed_facets[..],
                &mut undecidable_points,
            );

            // println!("Verifying facets at iteration: {}, k: {}", i, k);
            // for i in 0..triangles.len() {
            //     if triangles[i].valid {
            //         super::check_facet_links(i, &triangles[..]);
            //     }
            // }
        }

        i = i + 1;
    }

    let mut idx = Vec::new();

    for facet in triangles.iter() {
        if facet.valid {
            idx.push([
                facet.pts[0] as u32,
                facet.pts[1] as u32,
                facet.pts[2] as u32,
            ]);
        }
    }

    let mut points = points.to_vec();
    utils::remove_unused_points(&mut points, &mut idx[..]);

    assert!(points.len() != 0, "Internal error: empty output mesh.");
    // super::check_convex_hull(&points, &idx);

    (points, idx)
}

fn compute_silhouette(
    facet: usize,
    indirect_id: usize,
    point: usize,
    out_facets_and_idx: &mut Vec<(usize, usize)>,
    points: &[Point3<Real>],
    removed_facets: &mut Vec<usize>,
    triangles: &mut [TriangleFacet],
) {
    if triangles[facet].valid {
        if !triangles[facet].order_independent_can_be_seen_by_point(point, points) {
            // println!("triangles: {}, valid: true, keep: true", facet);
            // println!(
            //     "Taking edge: [{}, {}]",
            //     triangles[facet].second_point_from_edge(indirect_id),
            //     triangles[facet].first_point_from_edge(indirect_id)
            // );
            out_facets_and_idx.push((facet, indirect_id));
        } else {
            triangles[facet].valid = false; // The facet must be removed from the convex hull.
            removed_facets.push(facet);
            // println!("triangles: {}, valid: true, keep: false", facet);

            compute_silhouette(
                triangles[facet].adj[(indirect_id + 1) % 3],
                triangles[facet].indirect_adj_id[(indirect_id + 1) % 3],
                point,
                out_facets_and_idx,
                points,
                removed_facets,
                triangles,
            );

            compute_silhouette(
                triangles[facet].adj[(indirect_id + 2) % 3],
                triangles[facet].indirect_adj_id[(indirect_id + 2) % 3],
                point,
                out_facets_and_idx,
                points,
                removed_facets,
                triangles,
            );
        }
    } else {
        // println!("triangles: {}, valid: false, keep: false", facet);
    }
}

fn fix_silhouette_topology(
    points: &[Point3<Real>],
    out_facets_and_idx: &mut Vec<(usize, usize)>,
    removed_facets: &mut Vec<usize>,
    triangles: &mut [TriangleFacet],
) {
    // FIXME: don't allocate this everytime.
    let mut workspace = vec![0; points.len()];
    let mut needs_fixing = false;

    // NOTE: we wore with the second_point_from_edge instead
    // of the first one, because when we traverse the silhouette
    // we see the second edge point before the first.
    for (facet, adj_id) in &*out_facets_and_idx {
        let p = triangles[*facet].second_point_from_edge(*adj_id);
        workspace[p] += 1;

        if workspace[p] > 1 {
            needs_fixing = true;
        }
    }

    // We detected a topological problem, i.e., we have
    // multiple loops.
    if needs_fixing {
        // First, we need to know which loop is the one we
        // need to keep.
        let mut loop_start = 0;
        for (facet, adj_id) in &*out_facets_and_idx {
            let p1 = points[triangles[*facet].second_point_from_edge(*adj_id)];
            let p2 = points[triangles[*facet].first_point_from_edge(*adj_id)];
            let supp = indexed_support_point_nth(
                &(p2 - p1),
                points,
                out_facets_and_idx
                    .iter()
                    .map(|(f, ai)| triangles[*f].second_point_from_edge(*ai)),
            )
            .unwrap();
            let selected = &out_facets_and_idx[supp];
            if workspace[triangles[selected.0].second_point_from_edge(selected.1)] == 1 {
                // This is a valid point to start with.
                loop_start = supp;
                break;
            }
        }

        let mut removing = None;
        let old_facets_and_idx = std::mem::replace(out_facets_and_idx, Vec::new());

        for i in 0..old_facets_and_idx.len() {
            let facet_id = (loop_start + i) % old_facets_and_idx.len();
            let (facet, adj_id) = old_facets_and_idx[facet_id];

            match removing {
                Some(p) => {
                    let p1 = triangles[facet].second_point_from_edge(adj_id);
                    if p == p1 {
                        removing = None;
                    }
                }
                _ => {
                    let p1 = triangles[facet].second_point_from_edge(adj_id);
                    if workspace[p1] > 1 {
                        removing = Some(p1);
                    }
                }
            }

            if removing.is_some() {
                if triangles[facet].valid {
                    triangles[facet].valid = false;
                    removed_facets.push(facet);
                }
            } else {
                out_facets_and_idx.push((facet, adj_id));
            }

            // // Debug
            // {
            //     let p1 = triangles[facet].second_point_from_edge(adj_id);
            //     let p2 = triangles[facet].first_point_from_edge(adj_id);
            //     if removing.is_some() {
            //         print!("/{}, {}\\ ", p1, p2);
            //     } else {
            //         print!("[{}, {}] ", p1, p2);
            //     }
            // }
        }

        // println!("");
    }
}

fn attach_and_push_facets(
    silhouette_loop_facets_and_idx: &[(usize, usize)],
    point: usize,
    points: &[Point3<Real>],
    triangles: &mut Vec<TriangleFacet>,
    removed_facets: &[usize],
    undecidable: &mut Vec<usize>,
) {
    // The silhouette is built to be in CCW order.
    let mut new_facets = Vec::with_capacity(silhouette_loop_facets_and_idx.len());

    // Create new facets.
    let mut adj_facet: usize;
    let mut indirect_id: usize;

    for i in 0..silhouette_loop_facets_and_idx.len() {
        adj_facet = silhouette_loop_facets_and_idx[i].0;
        indirect_id = silhouette_loop_facets_and_idx[i].1;

        // print!(
        //     "[{}, {}] ",
        //     triangles[adj_facet].second_point_from_edge(indirect_id),
        //     triangles[adj_facet].first_point_from_edge(indirect_id)
        // );

        let facet = TriangleFacet::new(
            point,
            triangles[adj_facet].second_point_from_edge(indirect_id),
            triangles[adj_facet].first_point_from_edge(indirect_id),
            points,
        );
        new_facets.push(facet);
    }
    // println!("");

    // Link the facets together.
    for i in 0..silhouette_loop_facets_and_idx.len() {
        let prev_facet;

        if i == 0 {
            prev_facet = triangles.len() + silhouette_loop_facets_and_idx.len() - 1;
        } else {
            prev_facet = triangles.len() + i - 1;
        }

        let (middle_facet, middle_id) = silhouette_loop_facets_and_idx[i];
        let next_facet = triangles.len() + (i + 1) % silhouette_loop_facets_and_idx.len();

        new_facets[i].set_facets_adjascency(prev_facet, middle_facet, next_facet, 2, middle_id, 0);
        assert!(!triangles[triangles[middle_facet].adj[middle_id]].valid); // Check that we are not overwriting a valid link.
        triangles[middle_facet].adj[middle_id] = triangles.len() + i; // The future id of curr_facet.
        triangles[middle_facet].indirect_adj_id[middle_id] = 1;
    }

    // Assign to each facets some of the points which can see it.
    // FIXME: refactor this with the others.
    for curr_facet in removed_facets.iter() {
        for visible_point in triangles[*curr_facet].visible_points.iter() {
            if points[*visible_point] == points[point] {
                continue;
            }

            let mut furthest = usize::max_value();
            let mut furthest_dist = 0.0;

            for (i, curr_facet) in new_facets.iter_mut().enumerate() {
                if !curr_facet.affinely_dependent {
                    let distance = curr_facet.distance_to_point(*visible_point, points);

                    if distance > furthest_dist {
                        furthest = i;
                        furthest_dist = distance;
                    }
                }
            }

            if furthest != usize::max_value() {
                if new_facets[furthest].can_see_point(*visible_point, points) {
                    new_facets[furthest].add_visible_point(*visible_point, points);
                }
            }

            // If none of the facet can be seen from the point, it is implicitly
            // deleted because it won't be referenced by any facet.
        }
    }

    // Try to assign collinear points to one of the new facets.
    let mut i = 0;

    while i != undecidable.len() {
        let mut furthest = usize::max_value();
        let mut furthest_dist = 0.0;
        let undecidable_point = undecidable[i];

        for (j, curr_facet) in new_facets.iter_mut().enumerate() {
            if curr_facet.can_see_point(undecidable_point, points) {
                let distance = curr_facet.distance_to_point(undecidable_point, points);

                if distance > furthest_dist {
                    furthest = j;
                    furthest_dist = distance;
                }
            }
        }

        if furthest != usize::max_value() {
            new_facets[furthest].add_visible_point(undecidable_point, points);
            let _ = undecidable.swap_remove(i);
        } else {
            i = i + 1;
        }
    }

    // Push facets.
    // FIXME: can we avoid the tmp vector `new_facets` ?
    triangles.append(&mut new_facets);
}

#[cfg(test)]
mod test {
    use crate::transformation;
    #[cfg(feature = "dim2")]
    use na::Point2;

    #[cfg(feature = "dim2")]
    #[test]
    fn test_simple_convex_hull() {
        let points = [
            Point2::new(4.723881f32, 3.597233),
            Point2::new(3.333363, 3.429991),
            Point2::new(3.137215, 2.812263),
        ];

        let chull = transformation::convex_hull(points.as_slice());

        assert!(chull.coords.len() == 3);
    }

    #[cfg(feature = "dim3")]
    #[test]
    fn test_ball_convex_hull() {
        use crate::shape::Ball;

        // This triggered a failure to an affinely dependent facet.
        let (points, _) = Ball::new(0.4).to_trimesh(20, 20);
        let (vertices, _) = transformation::convex_hull(points.as_slice());

        // dummy test, we are just checking that the construction did not fail.
        assert!(vertices.len() == vertices.len());
    }
}
