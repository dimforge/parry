//! Covering a shape with the cells of a body-centered cubic lattice (the lattice of Labelle and
//! Shewchuk's isosurface stuffing, SIGGRAPH 2007, kept whole instead of cut). Every cell the shape
//! reaches is kept by exact intersection tests, so the mesh always contains what it covers.

use super::{MeshEnclosure, VolumeMesh, VolumeMeshParameters};
use crate::bounding_volume::{Aabb, BoundingVolume};
#[cfg(not(feature = "std"))]
use crate::math::ComplexField;
use crate::math::{Pose, Real, Vector};
use crate::query::PointQueryWithLocation;
use crate::shape::{SupportMap, TriMesh, TriMeshFlags};
use crate::utils::hashmap::HashMap;
use alloc::vec::Vec;

/// An order-preserving map over a slice, in parallel when the `parallel` feature is on.
///
/// Everything the cover pipeline hands this is a pure per-item computation, and rayon's
/// collect keeps item order, so the parallel and sequential paths produce identical
/// results, bit for bit.
#[cfg(feature = "parallel")]
pub(super) fn par_map<T: Sync, R: Send, F: Fn(&T) -> R + Sync + Send>(items: &[T], f: F) -> Vec<R> {
    use rayon::prelude::*;
    items.par_iter().map(f).collect()
}

/// The sequential twin of the parallel [`par_map`].
#[cfg(not(feature = "parallel"))]
pub(super) fn par_map<T, R, F: Fn(&T) -> R>(items: &[T], f: F) -> Vec<R> {
    items.iter().map(f).collect()
}

/// A convex point set as a support map, which is all GJK asks of a shape: the tetrahedra
/// and convex hulls the enclosure tests intersect never need more structure than this.
pub(super) struct PointsSupportMap<'a>(pub &'a [Vector]);

impl SupportMap for PointsSupportMap<'_> {
    fn local_support_point(&self, dir: Vector) -> Vector {
        let mut best = self.0[0];
        for pt in &self.0[1..] {
            if pt.dot(dir) > best.dot(dir) {
                best = *pt;
            }
        }
        best
    }
}

/// Whether the hull of a point set (a tetrahedron, a cage triangle) and a convex
/// support-mapped shape intersect (GJK).
pub(super) fn convex_points_intersect(points: &[Vector], other: &impl SupportMap) -> bool {
    crate::query::details::intersection_test_support_map_support_map(
        &Pose::IDENTITY,
        &PointsSupportMap(points),
        other,
    )
}

/// The oracle of [`tetrahedralize`]: the boundary mesh's own pseudo-normal signed distance
/// and exact intersection tests, which is everything the lattice reads of the shape.
pub(super) struct MeshOracle<'a> {
    mesh: &'a TriMesh,
}

impl<'a> MeshOracle<'a> {
    pub fn new(mesh: &'a TriMesh) -> Self {
        Self { mesh }
    }
}

impl MeshOracle<'_> {
    /// The signed distance to the boundary: positive inside.
    pub fn signed_distance(&self, pt: Vector) -> Real {
        signed_distance(self.mesh, pt)
    }

    /// Whether the shape's surface reaches into the box: exact (the BVH prunes to candidate
    /// triangles, GJK decides each, the first hit exits), which is what the cover's octree needs
    /// to decide a subdivision.
    pub fn crosses_region(&self, region: &Aabb) -> bool {
        let corners = [
            region.mins,
            Vector::new(region.maxs.x, region.mins.y, region.mins.z),
            Vector::new(region.mins.x, region.maxs.y, region.mins.z),
            Vector::new(region.maxs.x, region.maxs.y, region.mins.z),
            Vector::new(region.mins.x, region.mins.y, region.maxs.z),
            Vector::new(region.maxs.x, region.mins.y, region.maxs.z),
            Vector::new(region.mins.x, region.maxs.y, region.maxs.z),
            region.maxs,
        ];
        let crossed = self
            .mesh
            .bvh()
            .intersect_aabb(region)
            .any(|tri| convex_points_intersect(&corners, &self.mesh.triangle(tri)));
        crossed
    }

    /// Whether the shape's surface passes through the tetrahedron ([`MeshEnclosure::Crust`]):
    /// unlike [`Self::intersects`], a buried tetrahedron is not reached and no orientation is
    /// read, so the crust stays hollow and an open mesh can be crusted.
    pub fn surface_intersects(&self, tet: &[Vector; 4]) -> bool {
        let aabb = Aabb::from_points(tet.iter().copied());
        // Only the triangles the tetrahedron's Aabb can see are worth testing.
        let crossed = self
            .mesh
            .bvh()
            .intersect_aabb(&aabb)
            .any(|tri| convex_points_intersect(tet, &self.mesh.triangle(tri)));
        crossed
    }

    /// The triangles reaching into a region, collected once per boundary vertex or face for the
    /// smoothing's `_among` queries, so the per-move tests scan a few local triangles instead of
    /// re-walking the BVH.
    pub fn collect(&self, region: &Aabb, out: &mut Vec<u32>) {
        out.extend(self.mesh.bvh().intersect_aabb(region));
    }

    /// Whether the shape forbids a cage triangle passing there, tested against the listed
    /// triangles only: containment breaks exactly when the cage's boundary crosses the
    /// shape, so this is the smoothing's barrier.
    pub fn crosses_among(&self, triangle: &[Vector; 3], primitives: &[u32]) -> bool {
        // The Aabb reject is what keeps the listed neighborhood cheap: GJK only runs on
        // the one or two triangles overlapping the cage face's box.
        let aabb = Aabb::from_points(triangle.iter().copied());
        primitives.iter().any(|tri| {
            let candidate = self.mesh.triangle(*tri);
            let lows = candidate.a.min(candidate.b).min(candidate.c);
            let highs = candidate.a.max(candidate.b).max(candidate.c);
            lows.cmple(aabb.maxs).all()
                && highs.cmpge(aabb.mins).all()
                && convex_points_intersect(triangle, &candidate)
        })
    }

    /// The point of the shape's boundary closest to `pt` among the listed triangles, or
    /// `None` if the list is empty.
    pub fn project_among(&self, pt: Vector, primitives: &[u32]) -> Option<Vector> {
        use crate::query::PointQuery;

        // A triangle whose box is already farther than the best projection cannot beat it.
        let mut best = None;
        let mut best_dist = Real::MAX;
        for tri in primitives {
            let candidate = self.mesh.triangle(*tri);
            let lows = candidate.a.min(candidate.b).min(candidate.c);
            let highs = candidate.a.max(candidate.b).max(candidate.c);
            if (pt.clamp(lows, highs) - pt).length() >= best_dist {
                continue;
            }
            let projection = candidate.project_local_point(pt, true).point;
            let dist = (projection - pt).length();
            if dist < best_dist {
                best = Some(projection);
                best_dist = dist;
            }
        }
        best
    }
}

/// The body-centered cubic lattice covering an Aabb: the cube corners first, then the cube centers.
///
/// Coordinates are doubled so that both sublattices are integral: the lattice point with
/// coordinates `c` sits at `origin + c * cell_size / 2`, corners having only even coordinates and
/// centers only odd ones.
struct Lattice {
    origin: Vector,
    half_cell: Real,
    dims: [i32; 3],
    /// Number of even (resp. odd) coordinates along each axis.
    num_even: [i32; 3],
    num_odd: [i32; 3],
    /// Number of cube corners, i.e. the index of the first cube center.
    corners_len: usize,
}

impl Lattice {
    fn covering(aabb: Aabb, cell_size: Real) -> Self {
        let half_cell = cell_size * 0.5;
        // One cell of margin all around, so that every edge crossing the boundary has both of its
        // endpoints sampled.
        let origin = aabb.mins - Vector::splat(cell_size);
        let dims = core::array::from_fn(|i| {
            let extent = aabb.maxs[i] - aabb.mins[i] + cell_size * 2.0;
            ((extent / half_cell).ceil() as i32).max(2)
        });
        let num_even = core::array::from_fn(|i: usize| dims[i] / 2 + 1);
        let num_odd = core::array::from_fn(|i: usize| (dims[i] + 1) / 2);
        let corners_len = (num_even[0] * num_even[1] * num_even[2]) as usize;

        Self {
            origin,
            half_cell,
            dims,
            num_even,
            num_odd,
            corners_len,
        }
    }

    fn len(&self) -> usize {
        self.corners_len + (self.num_odd[0] * self.num_odd[1] * self.num_odd[2]) as usize
    }

    /// The index of the lattice point with the given doubled coordinates, if it is one and it is
    /// in range.
    fn index(&self, c: [i32; 3]) -> Option<u32> {
        if (0..3).any(|i| c[i] < 0 || c[i] > self.dims[i]) {
            return None;
        }

        let parity = c[0].rem_euclid(2);
        if c[1].rem_euclid(2) != parity || c[2].rem_euclid(2) != parity {
            return None;
        }

        let (base, counts) = if parity == 0 {
            (0, self.num_even)
        } else {
            (self.corners_len as i32, self.num_odd)
        };
        let [i, j, k] = [c[0] / 2, c[1] / 2, c[2] / 2];
        Some((base + i + (j + k * counts[1]) * counts[0]) as u32)
    }

    fn point(&self, c: [i32; 3]) -> Vector {
        self.origin + Vector::new(c[0] as Real, c[1] as Real, c[2] as Real) * self.half_cell
    }

    /// The doubled coordinates of every lattice point, in index order.
    fn coords(&self) -> Vec<[i32; 3]> {
        let mut coords = Vec::with_capacity(self.len());
        for (parity, counts) in [(0, self.num_even), (1, self.num_odd)] {
            for k in 0..counts[2] {
                for j in 0..counts[1] {
                    for i in 0..counts[0] {
                        coords.push([i * 2 + parity, j * 2 + parity, k * 2 + parity]);
                    }
                }
            }
        }
        coords
    }
}

/// The lattice tetrahedra, each of them stored so that its two long edges are `[0, 1]` and
/// `[2, 3]`.
///
/// The lattice decomposes into octahedra: two neighboring cube centers plus the four corners of the
/// cube face between them. Splitting each octahedron along its two centers gives four tetrahedra,
/// one per edge of that face.
fn lattice_tetrahedra(lattice: &Lattice, coords: &[[i32; 3]]) -> Vec<[u32; 4]> {
    let mut tets = Vec::new();

    for (center, &c) in coords.iter().enumerate().skip(lattice.corners_len) {
        let center = center as u32;

        for axis in 0..3 {
            let mut step = [0; 3];
            step[axis] = 2;
            let Some(opposite) = lattice.index(core::array::from_fn(|i| c[i] + step[i])) else {
                continue;
            };

            // The corners of the cube face between the two centers, in cyclic order.
            let (u, v) = ((axis + 1) % 3, (axis + 2) % 3);
            let mut corners = [0; 4];
            let mut complete = true;

            for (k, [su, sv]) in [[1, 1], [-1, 1], [-1, -1], [1, -1]].into_iter().enumerate() {
                let mut corner = c;
                corner[axis] += 1;
                corner[u] += su;
                corner[v] += sv;

                match lattice.index(corner) {
                    Some(corner) => corners[k] = corner,
                    None => {
                        complete = false;
                        break;
                    }
                }
            }

            if !complete {
                continue;
            }

            for k in 0..4 {
                tets.push([center, opposite, corners[k], corners[(k + 1) % 4]]);
            }
        }
    }

    tets
}

/// Whether every edge of the mesh is shared by exactly two triangles.
fn is_closed(indices: &[[u32; 3]]) -> bool {
    let mut edges: HashMap<[u32; 2], u32> = HashMap::default();

    for tri in indices {
        for k in 0..3 {
            *edges.entry(edge_key(tri[k], tri[(k + 1) % 3])).or_insert(0) += 1;
        }
    }

    edges.values().all(|count| *count == 2)
}

/// The signed distance to the boundary mesh: positive inside.
fn signed_distance(mesh: &TriMesh, pt: Vector) -> Real {
    // The projection has to keep its location for the inside test to use the pseudo-normals.
    let (proj, _) = mesh.project_local_point_and_get_location(pt, false);
    let dist = (pt - proj.point).length();
    if proj.is_inside {
        dist
    } else {
        -dist
    }
}

fn edge_key(a: u32, b: u32) -> [u32; 2] {
    if a < b {
        [a, b]
    } else {
        [b, a]
    }
}

/// The tetrahedron flipped positive, or nothing if it is degenerate.
fn orient_tet(pts: &[Vector], mut tet: [u32; 4], min_volume: Real) -> Option<[u32; 4]> {
    let [a, b, c, d] = tet.map(|i| pts[i as usize]);
    let volume = (b - a).cross(c - a).dot(d - a) / 6.0;

    if volume.abs() < min_volume {
        return None;
    }

    if volume < 0.0 {
        tet.swap(0, 1);
    }

    Some(tet)
}

/// A grid of cells filling the domain: the uniform lattice, or the cover's subdivided
/// octree grid.
pub(super) struct BackgroundGrid {
    pub points: Vec<Vector>,
    pub cells: Vec<[u32; 4]>,
}

/// The domain the grid has to cover: the shape's Aabb, grown by the whole cells the cover
/// keeps around the shape.
pub(super) fn domain(aabb: Aabb, cell_size: Real) -> Aabb {
    aabb.loosened(cell_size)
}

/// The body-centered cubic lattice covering an Aabb, as a background grid.
pub(super) fn uniform_grid(aabb: Aabb, cell_size: Real) -> BackgroundGrid {
    let lattice = Lattice::covering(aabb, cell_size);
    let coords = lattice.coords();
    let points = coords.iter().map(|&c| lattice.point(c)).collect();
    let cells = lattice_tetrahedra(&lattice, &coords);

    BackgroundGrid { points, cells }
}

/// Which non-crossing cells are inside the shape, by connectivity: cells untouched by the surface
/// that share a vertex lie on the same side of it, so their vertex-connected components are inside
/// or outside as a whole and one signed-distance probe per component classifies every cell in it.
fn flood_fill_inside(
    oracle: &MeshOracle,
    points: &[Vector],
    grid_cells: &[[u32; 4]],
    crossing: &[bool],
) -> Vec<bool> {
    // Union-find over the non-crossing cells, joined through shared vertices.
    let mut parent: Vec<u32> = (0..grid_cells.len() as u32).collect();
    fn root(parent: &mut [u32], mut i: u32) -> u32 {
        while parent[i as usize] != i {
            parent[i as usize] = parent[parent[i as usize] as usize];
            i = parent[i as usize];
        }
        i
    }

    let mut last_at_vertex = alloc::vec![u32::MAX; points.len()];
    for (id, cell) in grid_cells.iter().enumerate() {
        if crossing[id] {
            continue;
        }
        for v in cell {
            let previous = last_at_vertex[*v as usize];
            last_at_vertex[*v as usize] = id as u32;
            if previous != u32::MAX {
                let (a, b) = (root(&mut parent, previous), root(&mut parent, id as u32));
                parent[a as usize] = b;
            }
        }
    }

    // One probe per component, at the centroid of its representative cell.
    let mut inside_root: HashMap<u32, bool> = HashMap::default();
    let mut inside = alloc::vec![false; grid_cells.len()];
    for id in 0..grid_cells.len() {
        if crossing[id] {
            continue;
        }
        let component = root(&mut parent, id as u32);
        let is_inside = *inside_root.entry(component).or_insert_with(|| {
            let centroid = grid_cells[component as usize]
                .iter()
                .map(|v| points[*v as usize])
                .sum::<Vector>()
                / 4.0;
            oracle.signed_distance(centroid) > 0.0
        });
        inside[id] = is_inside;
    }

    inside
}

/// Keeps every whole cell of the background grid the shape reaches, cutting and warping
/// nothing: the mesh contains the shape by construction ([`MeshEnclosure::Cover`]).
///
/// A cell is kept when a vertex of it is inside the (dilated) shape, or when the exact
/// intersection test says the shape passes through it; the latter is what catches a feature
/// that slips between the sampled vertices, which the sampled field alone would lose.
fn cover_grid(
    grid: BackgroundGrid,
    oracle: &MeshOracle,
    params: &VolumeMeshParameters,
    min_volume: Real,
) -> Option<VolumeMesh> {
    let BackgroundGrid {
        points,
        cells: grid_cells,
    } = grid;

    // The cells the surface exactly passes through: the BVH prunes to the band along the
    // boundary, so the convex tests are only ever paid there.
    let crossing: Vec<bool> = par_map(&grid_cells, |tet| {
        let pts = tet.map(|v| points[v as usize]);
        oracle.surface_intersects(&pts)
    });

    // What else is kept: nothing for the crust (hollow by design); for the cover, the
    // cells inside the shape, classified by connectivity with one signed-distance probe
    // per component.
    let keep: Vec<bool> = if params.enclosure == MeshEnclosure::Crust {
        crossing
    } else {
        let inside = flood_fill_inside(oracle, &points, &grid_cells, &crossing);
        crossing
            .iter()
            .zip(&inside)
            .map(|(crossing, inside)| *crossing || *inside)
            .collect()
    };

    let ids: Vec<u32> = (0..grid_cells.len() as u32).collect();
    let kept = par_map(&ids, |id| {
        if keep[*id as usize] {
            orient_tet(&points, grid_cells[*id as usize], min_volume)
        } else {
            None
        }
    });
    let cells: Vec<[u32; 4]> = kept.into_iter().flatten().collect();

    if cells.is_empty() {
        return None;
    }

    let mut result = VolumeMesh {
        vertices: points,
        cells,
    };
    result.compact();

    Some(result)
}

pub fn tetrahedralize(
    vertices: &[Vector],
    indices: &[[u32; 3]],
    params: &VolumeMeshParameters,
) -> Option<VolumeMesh> {
    let cell_size = params.cell_size;

    if vertices.is_empty() || indices.is_empty() || cell_size <= 0.0 || cell_size.is_nan() {
        return None;
    }

    // The cut function is the signed distance to the boundary, so the boundary must be closed and
    // consistently oriented (duplicated vertices are welded first); the crust only keeps the cells
    // the surface passes through and never reads the sign, so it accepts an open mesh.
    let crust = params.enclosure == MeshEnclosure::Crust;
    let flags = if crust {
        TriMeshFlags::MERGE_DUPLICATE_VERTICES
    } else {
        TriMeshFlags::ORIENTED | TriMeshFlags::MERGE_DUPLICATE_VERTICES
    };
    let mesh = TriMesh::with_flags(vertices.to_vec(), indices.to_vec(), flags)
        .ok()
        .filter(|mesh| crust || (mesh.pseudo_normals().is_some() && is_closed(mesh.indices())))?;

    let oracle = MeshOracle::new(&mesh);
    let grid = if params.cover_subdivisions > 0 {
        super::cover_octree::cover_octree_grid(
            &oracle,
            domain(mesh.local_aabb(), cell_size),
            params,
        )?
    } else {
        uniform_grid(domain(mesh.local_aabb(), cell_size), cell_size)
    };

    let mut result = cover_grid(
        grid,
        &oracle,
        params,
        cell_size * cell_size * cell_size * 1.0e-6,
    )?;
    super::cover_smoothing::smooth_cover(&mut result, &oracle, params);
    Some(result)
}
