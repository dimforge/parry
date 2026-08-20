//! Triangulation of the interior of a closed polyline, by Delaunay refinement (Ruppert's and
//! Chew's algorithms, as implemented by `spade`).

use super::{VolumeMesh, VolumeMeshParameters};
use crate::math::{Real, Vector};
use crate::utils::hashmap::{Entry, HashMap};
use crate::utils::sanitize_spade_point;
use alloc::vec::Vec;
use spade::{
    AngleLimit, ConstrainedDelaunayTriangulation, Point2, RefinementParameters, Triangulation as _,
};

pub fn triangulate(
    vertices: &[Vector],
    indices: &[[u32; 2]],
    params: &VolumeMeshParameters,
) -> Option<VolumeMesh> {
    let cell_size = params.cell_size;

    if vertices.is_empty() || indices.is_empty() || cell_size <= 0.0 || cell_size.is_nan() {
        return None;
    }

    /*
     * The boundary, as constraint edges: only the vertices it references take part in the
     * triangulation.
     */
    let mut cdt: ConstrainedDelaunayTriangulation<Point2<Real>> =
        ConstrainedDelaunayTriangulation::new();
    let mut handles = HashMap::default();

    for edge in indices {
        let mut endpoints = [None; 2];

        for (k, vid) in edge.iter().enumerate() {
            let pt = *vertices.get(*vid as usize)?;
            let handle = match handles.entry(*vid) {
                Entry::Occupied(entry) => *entry.get(),
                Entry::Vacant(entry) => {
                    let pt = sanitize_spade_point(Point2::new(pt.x, pt.y));
                    *entry.insert(cdt.insert(pt).ok()?)
                }
            };
            endpoints[k] = Some(handle);
        }

        if let [Some(from), Some(to)] = endpoints {
            if from != to {
                // Intersecting constraint edges are left out instead of panicking; the boundary
                // shouldn't self-intersect in the first place.
                let _ = cdt.try_add_constraint(from, to);
            }
        }
    }

    /*
     * Refinement: split the triangles that are too large or too sharp, and the boundary edges that
     * stand in the way, until every triangle inside the boundary is fit for simulation.
     */
    let angle_limit = AngleLimit::from_rad(params.min_angle as f64);
    // The area of the equilateral triangle of side `cell_size`.
    let max_area = cell_size * cell_size * Real::sqrt(0.75) / 2.0;
    // The refinement gives up after a fixed number of added points, which defaults to ten times
    // the boundary's: far too few for a small `cell_size`, and a mesh that stops mid-refinement is
    // worse than a coarse one. Budget for a few times the elements the area target asks for.
    let enclosed_area = indices
        .iter()
        .map(|e| vertices[e[0] as usize].perp_dot(vertices[e[1] as usize]) / 2.0)
        .sum::<Real>()
        .abs();
    let budget = (enclosed_area / max_area * 4.0) as usize + 1000;
    let refinement = cdt.refine(
        RefinementParameters::new()
            .with_angle_limit(angle_limit)
            .with_max_allowed_area(max_area)
            .with_max_additional_vertices(budget)
            .exclude_outer_faces(true),
    );

    if !refinement.refinement_complete {
        log::warn!(
            "volume_mesh: the refinement ran out of vertices after {budget}; some elements are \
             badly shaped. Try a larger cell size or a smaller minimum angle."
        );
    }

    let outer: HashMap<_, ()> = refinement
        .excluded_faces
        .into_iter()
        .map(|face| (face.index(), ()))
        .collect();
    let mut vertices = alloc::vec![Vector::ZERO; cdt.num_vertices()];

    for vertex in cdt.vertices() {
        let pt = vertex.position();
        vertices[vertex.fix().index()] = Vector::new(pt.x, pt.y);
    }

    let mut cells = Vec::new();

    for face in cdt.inner_faces() {
        if outer.contains_key(&face.fix().index()) {
            continue;
        }

        let cell = face.vertices().map(|v| v.fix().index() as u32);
        let [a, b, c] = cell.map(|v| vertices[v as usize]);

        // Positively oriented, and no degenerate cell.
        let area = (b - a).perp_dot(c - a) / 2.0;
        if area.abs() < max_area * 1.0e-6 {
            continue;
        }

        cells.push(if area > 0.0 {
            cell
        } else {
            [cell[0], cell[2], cell[1]]
        });
    }

    if cells.is_empty() {
        return None;
    }

    let mut result = VolumeMesh { vertices, cells };
    result.compact();

    Some(result)
}
