//! Filling a closed boundary with a simulation-grade simplex mesh.

use crate::math::{Real, Vector, DIM};
use alloc::vec::Vec;

#[cfg(feature = "dim3")]
mod cover_octree;
#[cfg(feature = "dim3")]
mod cover_smoothing;
#[cfg(feature = "dim2")]
mod delaunay_refinement;
#[cfg(feature = "dim3")]
mod isosurface_stuffing;

/// A simplicial mesh filling the interior of a closed boundary: triangles in 2D, tetrahedra in 3D.
#[derive(Clone, Debug, Default)]
#[cfg_attr(
    feature = "serde-serialize",
    derive(serde::Serialize, serde::Deserialize)
)]
pub struct VolumeMesh {
    /// The mesh vertices.
    pub vertices: Vec<Vector>,
    /// The mesh cells, positively oriented (positive signed area/volume).
    pub cells: Vec<[u32; DIM + 1]>,
}

/// Which cover the mesh is (3D only): of the shape's volume, or of its surface alone.
///
/// Either way the mesh is made of whole lattice cells the shape reaches, so it *contains*
/// what it covers, which is what an embedding wants: every point of the shape (a skin
/// vertex, say) interpolates inside a cell instead of extrapolating outside all of them.
#[cfg(feature = "dim3")]
#[derive(Copy, Clone, Debug, PartialEq, Eq, Default)]
pub enum MeshEnclosure {
    /// Every lattice cell the shape reaches is kept whole: a solid fill containing the shape by
    /// construction, blocky at the cell size (see [`VolumeMeshParameters::cover_smoothing`] and
    /// [`VolumeMeshParameters::cover_subdivisions`]). Needs a closed, consistently oriented mesh.
    #[default]
    Cover,
    /// The cover of the shape's surface alone: only the cells the surface crosses are kept,
    /// the interior stays empty, and the mesh need not be closed or oriented. The result is a
    /// hollow shell that deforms like a shell, not like a solid.
    Crust,
}

/// Parameters of [`volume_mesh`].
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct VolumeMeshParameters {
    /// Target size of the generated elements.
    pub cell_size: Real,
    /// Smallest angle, in radians, the refinement aims for (2D only).
    ///
    /// Above 30 degrees the refinement is not guaranteed to terminate, so it may stop early and
    /// leave some elements badly shaped.
    #[cfg(feature = "dim2")]
    pub min_angle: Real,
    /// Whether the whole shape is covered or its surface alone (3D only); see
    /// [`MeshEnclosure`].
    #[cfg(feature = "dim3")]
    pub enclosure: MeshEnclosure,
    /// How many shrink-wrap iterations smooth the staircase of a [`MeshEnclosure::Cover`] mesh,
    /// zero leaving it raw (3D only, cover only). Each iteration pulls the boundary toward
    /// the shape, held off by [`Self::cover_guard`], and commits a move only if containment holds.
    #[cfg(feature = "dim3")]
    pub cover_smoothing: u32,
    /// How close to the shape the smoothed cover may pull its boundary, as a fraction of
    /// the local cell size (3D only, read by the cover's smoothing only).
    #[cfg(feature = "dim3")]
    pub cover_guard: Real,
    /// How many halvings below [`Self::cell_size`] a [`MeshEnclosure::Cover`] cell crossing the
    /// shape's boundary may be refined, zero keeping the boundary at the cell size (3D only, cover
    /// only). The refinement is on the octree (no hanging vertex) and runs before the smoothing.
    #[cfg(feature = "dim3")]
    pub cover_subdivisions: u32,
}

impl VolumeMeshParameters {
    /// Parameters generating elements of size `cell_size`: in 3D a raw cover (no
    /// smoothing, no subdivision), in 2D a refinement aiming for 30 degree angles.
    pub fn new(cell_size: Real) -> Self {
        Self {
            cell_size,
            #[cfg(feature = "dim2")]
            #[cfg_attr(feature = "f64", expect(clippy::unnecessary_cast))]
            min_angle: core::f64::consts::PI as Real / 6.0,
            #[cfg(feature = "dim3")]
            enclosure: MeshEnclosure::Cover,
            #[cfg(feature = "dim3")]
            cover_smoothing: 0,
            #[cfg(feature = "dim3")]
            cover_guard: 0.15,
            #[cfg(feature = "dim3")]
            cover_subdivisions: 0,
        }
    }
}

/// Fills a boundary `(vertices, indices)` (2D polyline, 3D triangle mesh) with a simplex mesh
/// for finite-element simulation: 2D refines a constrained Delaunay triangulation of the boundary,
/// 3D builds a lattice cover containing it. Returns `None` if the boundary encloses nothing.
pub fn volume_mesh(
    vertices: &[Vector],
    indices: &[[u32; DIM]],
    params: &VolumeMeshParameters,
) -> Option<VolumeMesh> {
    #[cfg(feature = "dim2")]
    {
        delaunay_refinement::triangulate(vertices, indices, params)
    }
    #[cfg(feature = "dim3")]
    {
        isosurface_stuffing::tetrahedralize(vertices, indices, params)
    }
}

impl VolumeMesh {
    /// The connected component of each cell: two cells sharing a vertex are in the same one.
    ///
    /// Components are numbered from zero, so the number of them is the largest index plus one.
    /// A mesh in more than one piece means the input's own surface is in several pieces: the
    /// cover keeps every cell the surface reaches, however thin the feature.
    pub fn connected_components(&self) -> Vec<u32> {
        let mut parent: Vec<u32> = (0..self.vertices.len() as u32).collect();

        fn root(parent: &mut [u32], mut i: u32) -> u32 {
            while parent[i as usize] != i {
                parent[i as usize] = parent[parent[i as usize] as usize];
                i = parent[i as usize];
            }
            i
        }

        for cell in &self.cells {
            for vid in &cell[1..] {
                let (a, b) = (root(&mut parent, cell[0]), root(&mut parent, *vid));
                parent[a as usize] = b;
            }
        }

        let mut ids = alloc::vec![u32::MAX; self.vertices.len()];
        let mut count = 0;

        self.cells
            .iter()
            .map(|cell| {
                let id = &mut ids[root(&mut parent, cell[0]) as usize];
                if *id == u32::MAX {
                    *id = count;
                    count += 1;
                }
                *id
            })
            .collect()
    }

    /// Removes the vertices no cell references, and reindexes the cells accordingly.
    pub fn compact(&mut self) {
        let mut remap = alloc::vec![u32::MAX; self.vertices.len()];
        let mut vertices = Vec::with_capacity(self.vertices.len());

        for cell in &mut self.cells {
            for vid in cell {
                let remapped = &mut remap[*vid as usize];
                if *remapped == u32::MAX {
                    *remapped = vertices.len() as u32;
                    vertices.push(self.vertices[*vid as usize]);
                }
                *vid = *remapped;
            }
        }

        self.vertices = vertices;
    }
}
