//! Shrink-wrap of a [`super::MeshEnclosure::Cover`] mesh: flattens its staircase while keeping its
//! containment guarantee, a boundary vertex moving only if its cells keep a volume floor and its
//! boundary faces still do not cross the shape. Vertices relax color by color in a fixed order.

use super::isosurface_stuffing::{par_map, MeshOracle};
use super::{VolumeMesh, VolumeMeshParameters};
use crate::bounding_volume::{Aabb, BoundingVolume};
#[cfg(not(feature = "std"))]
use crate::math::ComplexField;
use crate::math::{Real, Vector};
use alloc::vec::Vec;

/// How much of its lattice volume a cell must keep for a move to commit.
const VOLUME_FLOOR: Real = 0.2;

/// How far a boundary vertex may travel from where the lattice put it, in local cell
/// sizes. The wrap's targets all lie within a cell, so this hardly ever binds; what it
/// buys is a bound that makes each vertex's neighborhood of the shape collectable once.
const TRAVEL_BUDGET: Real = 1.0;

/// The half-extent of the region a vertex's shape primitives are collected in, in local
/// cell sizes: the incident faces' own extent, plus every incident vertex's travel budget,
/// with margin.
const REGION: Real = 2.5;

/// The committed step, as a fraction of the guard, under which a vertex goes to sleep.
const SLEEP: Real = 0.05;

/// The four faces of a tetrahedron, as index triples into the cell.
const FACES: [[usize; 3]; 4] = [[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]];

fn cell_volume(positions: &[Vector], cell: [u32; 4]) -> Real {
    let [a, b, c, d] = cell.map(|i| positions[i as usize]);
    (b - a).cross(c - a).dot(d - a) / 6.0
}

/// The cell's volume, with vertex `v` read at `at` instead of its stored position.
fn cell_volume_at(positions: &[Vector], cell: [u32; 4], v: u32, at: Vector) -> Real {
    let p = cell.map(|i| if i == v { at } else { positions[i as usize] });
    (p[1] - p[0]).cross(p[2] - p[0]).dot(p[3] - p[0]) / 6.0
}

/// The face's points, with vertex `v` read at `at` instead of its stored position.
fn face_at(positions: &[Vector], face: [u32; 3], v: u32, at: Vector) -> [Vector; 3] {
    face.map(|i| if i == v { at } else { positions[i as usize] })
}

pub(super) fn smooth_cover(
    mesh: &mut VolumeMesh,
    oracle: &MeshOracle,
    params: &VolumeMeshParameters,
) {
    if params.cover_smoothing == 0 || params.cover_guard <= 0.0 || mesh.cells.is_empty() {
        return;
    }

    /*
     * The cover's boundary: the faces one cell shares with no other (found by sorting the
     * face keys, which is parallel-friendly and gives a deterministic order), their
     * vertices, and the adjacency the relaxation reads.
     */
    let mut keys: Vec<[u32; 3]> = Vec::with_capacity(mesh.cells.len() * 4);
    for cell in &mesh.cells {
        for face in FACES {
            let mut key = face.map(|k| cell[k]);
            key.sort_unstable();
            keys.push(key);
        }
    }
    #[cfg(feature = "parallel")]
    {
        use rayon::prelude::*;
        keys.par_sort_unstable();
    }
    #[cfg(not(feature = "parallel"))]
    keys.sort_unstable();

    let mut boundary_faces: Vec<[u32; 3]> = Vec::new();
    let mut at = 0;
    while at < keys.len() {
        let mut next = at + 1;
        while next < keys.len() && keys[next] == keys[at] {
            next += 1;
        }
        if next - at == 1 {
            boundary_faces.push(keys[at]);
        }
        at = next;
    }

    let n = mesh.vertices.len();
    let mut is_boundary = alloc::vec![false; n];
    let mut vertex_faces: Vec<Vec<u32>> = alloc::vec![Vec::new(); n];
    let mut neighbors: Vec<Vec<u32>> = alloc::vec![Vec::new(); n];

    for (id, face) in boundary_faces.iter().enumerate() {
        for k in 0..3 {
            let v = face[k] as usize;
            is_boundary[v] = true;
            vertex_faces[v].push(id as u32);
            for other in [face[(k + 1) % 3], face[(k + 2) % 3]] {
                if !neighbors[v].contains(&other) {
                    neighbors[v].push(other);
                }
            }
        }
    }

    let mut vertex_cells: Vec<Vec<u32>> = alloc::vec![Vec::new(); n];
    for (id, cell) in mesh.cells.iter().enumerate() {
        for v in cell {
            vertex_cells[*v as usize].push(id as u32);
        }
    }

    /*
     * What the pristine lattice grants each cell and vertex: the volume floor, and the
     * guard and step cap scaled by the local cell size.
     */
    let floors: Vec<Real> = mesh
        .cells
        .iter()
        .map(|cell| cell_volume(&mesh.vertices, *cell) * VOLUME_FLOOR)
        .collect();
    let scales: Vec<Real> = (0..n)
        .map(|v| {
            if !is_boundary[v] {
                return 0.0;
            }
            vertex_cells[v]
                .iter()
                .flat_map(|c| {
                    let pts = mesh.cells[*c as usize].map(|i| mesh.vertices[i as usize]);
                    (0..4).flat_map(move |i| (i + 1..4).map(move |j| (pts[i] - pts[j]).length()))
                })
                .fold(0.0, Real::max)
        })
        .collect();
    let guards: Vec<Real> = scales.iter().map(|s| params.cover_guard * s).collect();

    let boundary_vertices: Vec<u32> = (0..n as u32).filter(|v| is_boundary[*v as usize]).collect();

    /*
     * Each vertex's and face's neighborhood of the shape, collected once: the shape is static and
     * vertex travel is bounded, so the per-move barrier and projection scan these short lists
     * instead of the whole shape.
     */
    let origins: Vec<Vector> = mesh.vertices.clone();
    let neighborhoods: Vec<Vec<u32>> = par_map(&boundary_vertices, |&v| {
        let vid = v as usize;
        let region = Aabb::from_half_extents(origins[vid], Vector::splat(scales[vid] * REGION));
        let mut primitives = Vec::new();
        oracle.collect(&region, &mut primitives);
        primitives
    });
    let neighborhood_of = {
        let mut ids = alloc::vec![u32::MAX; n];
        for (k, v) in boundary_vertices.iter().enumerate() {
            ids[*v as usize] = k as u32;
        }
        ids
    };
    let face_lists: Vec<Vec<u32>> = par_map(&boundary_faces, |face| {
        let sweep = face
            .iter()
            .map(|v| scales[*v as usize])
            .fold(0.0, Real::max)
            * (TRAVEL_BUDGET * 1.1);
        let region = Aabb::from_points(face.iter().map(|v| origins[*v as usize])).loosened(sweep);
        let mut primitives = Vec::new();
        oracle.collect(&region, &mut primitives);
        primitives
    });

    /*
     * Greedy coloring of boundary vertices, no two of a color sharing a cell: a vertex's target
     * and validity read only its incident cells and faces, so each color is one parallel batch.
     */
    let mut colors = alloc::vec![u32::MAX; n];
    let mut palette: Vec<Vec<u32>> = Vec::new();
    for &v in &boundary_vertices {
        let mut used = alloc::vec![false; palette.len()];
        for &c in &vertex_cells[v as usize] {
            for o in mesh.cells[c as usize] {
                let color = colors[o as usize];
                if color != u32::MAX {
                    used[color as usize] = true;
                }
            }
        }
        let color = used.iter().position(|u| !*u).unwrap_or(palette.len());
        if color == palette.len() {
            palette.push(Vec::new());
        }
        colors[v as usize] = color as u32;
        palette[color].push(v);
    }

    /*
     * The wrap: color by color, each vertex validates its move against the batch's starting state
     * and proposals apply in order; the step is capped at the guard and tested at its midpoint, so
     * no sweep exceeds half a guard between exact tests; idle vertices wake on a neighbor's move.
     */
    let mut active = is_boundary.clone();

    for _ in 0..params.cover_smoothing {
        let mut moved_any = false;

        for group in &palette {
            let members: Vec<u32> = group
                .iter()
                .copied()
                .filter(|v| active[*v as usize])
                .collect();
            if members.is_empty() {
                continue;
            }

            let positions = &mesh.vertices;
            let proposals: Vec<(Vector, Real)> = par_map(&members, |&v| {
                let vid = v as usize;
                let pos = positions[vid];
                let guard = guards[vid];
                let primitives = &neighborhoods[neighborhood_of[vid] as usize];
                if neighbors[vid].is_empty() || primitives.is_empty() {
                    return (pos, 0.0);
                }

                let average = neighbors[vid]
                    .iter()
                    .map(|o| positions[*o as usize])
                    .sum::<Vector>()
                    / neighbors[vid].len() as Real;
                let projection = oracle.project_among(pos, primitives).unwrap_or(pos);
                let Some(outward) = (pos - projection).try_normalize() else {
                    return (pos, 0.0);
                };
                let held_off = projection + outward * guard;
                let target = (average + held_off) * 0.5;

                // The mid-step barrier test keeps the effective sweep at half of this, so
                // a full-guard cap converges twice as fast at the same tunneling
                // granularity.
                let mut step = target - pos;
                if step.length() > guard {
                    step = step.normalize() * guard;
                }
                // The travel budget: what makes the collected neighborhoods complete.
                let budget = scales[vid] * TRAVEL_BUDGET;
                let strayed = pos + step - origins[vid];
                if strayed.length() > budget {
                    step = origins[vid] + strayed.normalize() * budget - pos;
                }
                // A step already below the sleep threshold has nothing to buy: skip the
                // barrier work it would cost.
                if step.length() < guard * SLEEP {
                    return (pos, 0.0);
                }

                let valid = |to: Vector| -> bool {
                    for &c in &vertex_cells[vid] {
                        if cell_volume_at(positions, mesh.cells[c as usize], v, to)
                            < floors[c as usize]
                        {
                            return false;
                        }
                    }
                    let middle = (pos + to) * 0.5;
                    for &f in &vertex_faces[vid] {
                        let face = boundary_faces[f as usize];
                        let barrier = &face_lists[f as usize];
                        if oracle.crosses_among(&face_at(positions, face, v, to), barrier)
                            || oracle.crosses_among(&face_at(positions, face, v, middle), barrier)
                        {
                            return false;
                        }
                    }
                    true
                };

                for factor in [1.0, 0.5, 0.25] {
                    let candidate = pos + step * factor;
                    if valid(candidate) {
                        return (candidate, step.length() * factor);
                    }
                }
                (pos, 0.0)
            });

            for (v, (to, committed)) in members.iter().zip(&proposals) {
                let vid = *v as usize;
                mesh.vertices[vid] = *to;
                if *committed < guards[vid] * SLEEP {
                    active[vid] = false;
                } else {
                    moved_any = true;
                    for o in &neighbors[vid] {
                        active[*o as usize] = true;
                    }
                }
            }
        }

        if !moved_any {
            break;
        }
    }

    /*
     * A light interior relaxation: the ring of interior vertices next to the wrapped boundary
     * absorbed the squeeze, so they are eased toward their neighbors' average; the volume floor
     * still guards the cage's region against these moves.
     */
    let ring: Vec<u32> = (0..n as u32)
        .filter(|v| {
            !is_boundary[*v as usize]
                && vertex_cells[*v as usize].iter().any(|c| {
                    mesh.cells[*c as usize]
                        .iter()
                        .any(|o| is_boundary[*o as usize])
                })
        })
        .collect();

    for _ in 0..2 {
        for &v in &ring {
            let vid = v as usize;
            let mut sum = Vector::ZERO;
            let mut count = 0;
            for &c in &vertex_cells[vid] {
                for o in mesh.cells[c as usize] {
                    if o != v {
                        sum += mesh.vertices[o as usize];
                        count += 1;
                    }
                }
            }
            if count == 0 {
                continue;
            }

            let pos = mesh.vertices[vid];
            let step = sum / count as Real - pos;
            for factor in [0.5, 0.25] {
                let candidate = pos + step * factor;
                mesh.vertices[vid] = candidate;
                if vertex_cells[vid].iter().all(|c| {
                    cell_volume(&mesh.vertices, mesh.cells[*c as usize]) >= floors[*c as usize]
                }) {
                    break;
                }
                mesh.vertices[vid] = pos;
            }
        }
    }
}
