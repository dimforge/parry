//! The cover's subdivided background grid: an octree refined below the cell size where octants
//! cross the shape's boundary, turned into a conforming grid of cells (no hanging vertex) by the
//! balance rule and the transition cells.

use super::isosurface_stuffing::{par_map, BackgroundGrid, MeshOracle};
use super::VolumeMeshParameters;
use crate::bounding_volume::Aabb;
#[cfg(not(feature = "std"))]
use crate::math::ComplexField;
use crate::math::{Real, Vector};
use crate::utils::hashmap::HashMap;
use alloc::vec::Vec;

/// The eight corners of an octant, as offsets in units of its width.
const CORNERS: [[i32; 3]; 8] = [
    [0, 0, 0],
    [1, 0, 0],
    [0, 1, 0],
    [1, 1, 0],
    [0, 0, 1],
    [1, 0, 1],
    [0, 1, 1],
    [1, 1, 1],
];

/// The octants the boundary passes through, refined down to the finest level, and the coarser ones
/// filling the rest.
///
/// Octant coordinates are indices at their own level; vertex coordinates are integers in units of
/// half the finest octant's width, which makes every corner, center, face center and edge midpoint
/// of every level integral.
struct Octree {
    origin: Vector,
    /// Half the width of a finest octant, the unit of the vertex coordinates.
    half: Real,
    /// The level of the coarsest octants: they are `1 << levels` finest octants wide.
    levels: u32,
    /// The number of coarsest octants along each axis.
    dims: [i32; 3],
    /// The octants that have children.
    subdivided: HashMap<(u32, [i32; 3]), ()>,
}

impl Octree {
    /// The width of an octant, in vertex coordinates: a finest octant is two of them wide.
    fn width(&self, level: u32) -> i32 {
        2 << level
    }

    /// The number of octants of a level along each axis.
    fn count(&self, level: u32) -> [i32; 3] {
        core::array::from_fn(|k| self.dims[k] << (self.levels - level))
    }

    fn in_range(&self, level: u32, c: [i32; 3]) -> bool {
        let count = self.count(level);
        (0..3).all(|k| c[k] >= 0 && c[k] < count[k])
    }

    fn is_subdivided(&self, level: u32, c: [i32; 3]) -> bool {
        self.subdivided.contains_key(&(level, c))
    }

    /// A corner of an octant, in vertex coordinates.
    fn corner(&self, level: u32, c: [i32; 3], offset: [i32; 3]) -> [i32; 3] {
        let width = self.width(level);
        core::array::from_fn(|k| (c[k] + offset[k]) * width)
    }

    /// The center of an octant, in vertex coordinates.
    fn center(&self, level: u32, c: [i32; 3]) -> [i32; 3] {
        let width = self.width(level);
        core::array::from_fn(|k| c[k] * width + width / 2)
    }

    fn point(&self, v: [i32; 3]) -> Vector {
        self.origin + Vector::new(v[0] as Real, v[1] as Real, v[2] as Real) * self.half
    }

    /// The leaf covering a given octant: itself when it is one, an ancestor when it sits inside a
    /// coarser leaf, `None` when it is subdivided (several leaves cover it) or out of range.
    fn covering_leaf(&self, level: u32, c: [i32; 3]) -> Option<(u32, [i32; 3])> {
        if !self.in_range(level, c) {
            return None;
        }

        let mut current = self.levels;
        loop {
            let shifted: [i32; 3] = core::array::from_fn(|k| c[k] >> (current - level));
            if !self.is_subdivided(current, shifted) {
                return Some((current, shifted));
            }
            if current == level {
                return None;
            }
            current -= 1;
        }
    }

    /// Gives an octant children, so that its own children become leaves.
    fn split(&mut self, level: u32, c: [i32; 3]) {
        debug_assert!(level > 0, "a finest octant has no children");
        let _ = self.subdivided.insert((level, c), ());
    }

    /// Every leaf of the octree.
    fn leaves(&self) -> Vec<(u32, [i32; 3])> {
        let mut leaves = Vec::new();
        let mut stack: Vec<(u32, [i32; 3])> = Vec::new();
        let roots = self.count(self.levels);

        for k in 0..roots[2] {
            for j in 0..roots[1] {
                for i in 0..roots[0] {
                    stack.push((self.levels, [i, j, k]));
                }
            }
        }

        while let Some((level, c)) = stack.pop() {
            if level > 0 && self.is_subdivided(level, c) {
                for offset in CORNERS {
                    stack.push((level - 1, core::array::from_fn(|k| c[k] * 2 + offset[k])));
                }
            } else {
                leaves.push((level, c));
            }
        }

        leaves
    }

    /// Whether an octant finer than `level` puts a vertex at the given coordinates.
    fn has_finer_vertex(&self, level: u32, v: [i32; 3]) -> bool {
        if level == 0 {
            return false;
        }

        let finer = level - 1;
        let width = self.width(finer);
        if (0..3).any(|k| v[k].rem_euclid(width) != 0) {
            return false;
        }

        let base: [i32; 3] = core::array::from_fn(|k| v[k] / width);
        CORNERS.iter().any(|offset| {
            let c: [i32; 3] = core::array::from_fn(|k| base[k] - offset[k]);
            matches!(self.covering_leaf(finer, c), Some((leaf, _)) if leaf <= finer)
        })
    }
}

/// The center of a face of an octant, in vertex coordinates.
fn face_center(octree: &Octree, level: u32, c: [i32; 3], axis: usize, positive: bool) -> [i32; 3] {
    let corners = face_corners(axis, positive).map(|offset| octree.corner(level, c, offset));
    core::array::from_fn(|k| (corners[0][k] + corners[2][k]) / 2)
}

/// The four corners of a face of an octant, walking around it, as offsets in units of its width.
fn face_corners(axis: usize, positive: bool) -> [[i32; 3]; 4] {
    let (u, v) = ((axis + 1) % 3, (axis + 2) % 3);
    let base = i32::from(positive);
    [[0, 0], [1, 0], [1, 1], [0, 1]].map(|[du, dv]| {
        let mut corner = [0; 3];
        corner[axis] = base;
        corner[u] = du;
        corner[v] = dv;
        corner
    })
}

/// The background grid of a cover ([`MeshEnclosure::Cover`](super::MeshEnclosure::Cover)): an
/// octant that may cross the shape's boundary is refined
/// [`VolumeMeshParameters::cover_subdivisions`] halvings below the cell size and one that cannot
/// is not; balancing and transition cells keep the grid conforming.
pub(super) fn cover_octree_grid(
    oracle: &MeshOracle,
    aabb: Aabb,
    params: &VolumeMeshParameters,
) -> Option<BackgroundGrid> {
    let cell_size = params.cell_size;
    let fine = cell_size / (1 << params.cover_subdivisions) as Real;
    let levels = params.cover_subdivisions;

    // The domain, in cell-sized octants, with margin so the boundary never reaches its
    // border.
    let coarse_width = cell_size;
    let margin = cell_size * 2.0;
    let origin = aabb.mins - Vector::splat(margin);
    let dims: [i32; 3] = core::array::from_fn(|k| {
        let extent = aabb.maxs[k] - aabb.mins[k] + margin * 2.0;
        ((extent / coarse_width).ceil() as i32).max(1)
    });

    let mut octree = Octree {
        origin,
        half: fine * 0.5,
        levels,
        dims,
        subdivided: HashMap::default(),
    };

    let roots = octree.count(levels);
    let mut frontier: Vec<[i32; 3]> = Vec::new();
    for k in 0..roots[2] {
        for j in 0..roots[1] {
            for i in 0..roots[0] {
                frontier.push([i, j, k]);
            }
        }
    }

    // Level-synchronous refinement: every octant of a level tests independently whether the
    // boundary reaches into it (an exact, early-exit BVH existence query), the splits are applied
    // in one pass, and the children become the next level's frontier.
    for level in (1..=levels).rev() {
        let decisions: Vec<bool> = par_map(&frontier, |c| {
            let octant = Aabb::new(
                octree.point(octree.corner(level, *c, [0, 0, 0])),
                octree.point(octree.corner(level, *c, [1, 1, 1])),
            );
            oracle.crosses_region(&octant)
        });

        let mut next = Vec::new();
        for (c, split) in frontier.iter().zip(&decisions) {
            if *split {
                octree.split(level, *c);
                for offset in CORNERS {
                    next.push(core::array::from_fn(|k| c[k] * 2 + offset[k]));
                }
            }
        }
        frontier = next;
        if frontier.is_empty() {
            break;
        }
    }

    balance(&mut octree);
    Some(background_grid(&octree))
}

/// The Weak Balance Condition: two octants sharing so much as an edge may not differ by more than
/// one level, or there is no bridging them.
fn balance(octree: &mut Octree) {
    loop {
        let leaves = octree.leaves();

        // Each round decides every leaf against the same snapshot of the octree, in
        // parallel; balancing is monotone (splitting only ever demands more splits), so
        // the rounds converge on the same closure the one-by-one sweep would.
        let decisions: Vec<bool> = par_map(&leaves, |&(level, c)| {
            if level == 0 {
                return false;
            }

            // Every finest octant one step outside this leaf must be covered by a leaf no more
            // than one level coarser than it. Only the shell around the leaf is walked: the
            // octants inside it are its own.
            let span = 1 << level;
            let mut split = false;

            'shell: for axis in 0..3 {
                for side in [-1, span] {
                    for a in -1..=span {
                        for b in -1..=span {
                            let mut offset = [0; 3];
                            offset[axis] = side;
                            offset[(axis + 1) % 3] = a;
                            offset[(axis + 2) % 3] = b;
                            let neighbor: [i32; 3] =
                                core::array::from_fn(|k| c[k] * span + offset[k]);

                            let Some((neighbor_level, _)) = octree.covering_leaf(0, neighbor)
                            else {
                                continue;
                            };
                            // A neighbor two or more levels finer: this leaf is the one that has
                            // to give, or nothing can bridge them.
                            if neighbor_level + 1 < level {
                                split = true;
                                break 'shell;
                            }
                        }
                    }
                }
            }

            split
        });

        let mut changed = false;
        for (&(level, c), split) in leaves.iter().zip(&decisions) {
            if *split {
                octree.split(level, c);
                changed = true;
            }
        }

        if !changed {
            break;
        }
    }
}

/// The background grid of a balanced octree: figure 10 of the paper.
///
/// A cell spanning two octants would be created from both of them; the same-size case, the only
/// symmetric one, is built from the positive side alone.
fn background_grid(octree: &Octree) -> BackgroundGrid {
    // Every leaf emits its cells independently, as vertex coordinates: the octree is only
    // read, so the leaves are one parallel pass, and the coordinates are deduplicated into
    // ids afterward.
    let leaves = octree.leaves();
    let per_leaf: Vec<Vec<[[i32; 3]; 4]>> =
        par_map(&leaves, |&(level, c)| leaf_cells(octree, level, c));
    let coord_cells: Vec<[[i32; 3]; 4]> = per_leaf.into_iter().flatten().collect();

    let mut coords: Vec<[i32; 3]> = coord_cells.iter().flatten().copied().collect();
    #[cfg(feature = "parallel")]
    {
        use rayon::prelude::*;
        coords.par_sort_unstable();
    }
    #[cfg(not(feature = "parallel"))]
    coords.sort_unstable();
    coords.dedup();

    let mut ids: HashMap<[i32; 3], u32> = HashMap::default();
    for (id, v) in coords.iter().enumerate() {
        let _ = ids.insert(*v, id as u32);
    }

    let points = par_map(&coords, |v| octree.point(*v));
    let cells = par_map(&coord_cells, |quad| quad.map(|v| ids[&v]));

    BackgroundGrid { points, cells }
}

/// The cells of one leaf, as vertex coordinates: figure 10 of Labelle and Shewchuk.
///
/// A cell spanning two octants would be created from both of them; the same-size case, the
/// only symmetric one, is built from the positive side alone.
fn leaf_cells(octree: &Octree, level: u32, c: [i32; 3]) -> Vec<[[i32; 3]; 4]> {
    let mut cells = Vec::new();
    let center = octree.center(level, c);

    for axis in 0..3 {
        for positive in [false, true] {
            let face: [[i32; 3]; 4] =
                face_corners(axis, positive).map(|offset| octree.corner(level, c, offset));

            let mut neighbor = c;
            neighbor[axis] += if positive { 1 } else { -1 };

            if octree.in_range(level, neighbor) && octree.is_subdivided(level, neighbor) {
                // Finer octants across the face put a vertex at its center: quadrisected cells.
                let middle: [i32; 3] = core::array::from_fn(|k| (face[0][k] + face[2][k]) / 2);

                for e in 0..4 {
                    let (a, b) = (face[e], face[(e + 1) % 4]);
                    let midpoint: [i32; 3] = core::array::from_fn(|k| (a[k] + b[k]) / 2);

                    if octree.has_finer_vertex(level, midpoint) {
                        cells.push([middle, center, a, midpoint]);
                        cells.push([middle, center, midpoint, b]);
                    } else {
                        cells.push([middle, center, a, b]);
                    }
                }
                continue;
            }

            match octree.covering_leaf(level, neighbor) {
                Some((neighbor_level, neighbor_c)) if neighbor_level == level => {
                    // The same size: the lattice cells, built once per shared face.
                    if !positive {
                        continue;
                    }
                    let opposite = octree.center(neighbor_level, neighbor_c);

                    for e in 0..4 {
                        let (a, b) = (face[e], face[(e + 1) % 4]);
                        let midpoint: [i32; 3] = core::array::from_fn(|k| (a[k] + b[k]) / 2);

                        if octree.has_finer_vertex(level, midpoint) {
                            // A finer octant split this edge: two bisected cells, which the
                            // boundary never crosses.
                            cells.push([center, opposite, a, midpoint]);
                            cells.push([center, opposite, midpoint, b]);
                        } else {
                            cells.push([center, opposite, a, b]);
                        }
                    }
                }
                // A coarser neighbor, or the domain's border: two half-pyramids over the face,
                // whose diagonal runs through the coarse face's center so that its two triangles
                // are among the eight the coarse side fans out from that center.
                other => {
                    let diagonal = other
                        .and_then(|(neighbor_level, neighbor_c)| {
                            let middle =
                                face_center(octree, neighbor_level, neighbor_c, axis, !positive);
                            face.iter().position(|v| *v == middle)
                        })
                        .unwrap_or(0);
                    let (a, b, cc, d) = (
                        face[diagonal],
                        face[(diagonal + 1) % 4],
                        face[(diagonal + 2) % 4],
                        face[(diagonal + 3) % 4],
                    );
                    cells.push([center, a, b, cc]);
                    cells.push([center, a, cc, d]);
                }
            }
        }
    }

    cells
}
