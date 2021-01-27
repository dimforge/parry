// Rust port, with modifications, of https://github.com/kmammou/v-hacd/blob/master/src/VHACD_Lib/src/vhacdVolume.cpp
// By Khaled Mamou
//
// # License of the original C++ code:
// > Copyright (c) 2011 Khaled Mamou (kmamou at gmail dot com)
// > All rights reserved.
// >
// >
// > Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
// >
// > 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
// >
// > 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
// >
// > 3. The names of the contributors may not be used to endorse or promote products derived from this software without specific prior written permission.
// >
// > THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

use super::{FillMode, VoxelizedVolume};
use crate::bounding_volume::AABB;
use crate::math::{Matrix, Point, Real, Vector, DIM};
use crate::transformation::vhacd::CutPlane;
use std::sync::Arc;

#[cfg(feature = "dim2")]
type ConvexHull = Vec<Point<Real>>;
#[cfg(feature = "dim3")]
type ConvexHull = (Vec<Point<Real>>, Vec<[u32; DIM]>);

/// A voxel.
#[derive(Copy, Clone, Debug)]
pub struct Voxel {
    /// The integer coordinates of the voxel as part of the voxel grid.
    pub coords: Point<u32>,
    /// Is this voxel on the surface of the volume (i.e. not inside of it)?
    pub is_on_surface: bool,
    /// Range of indices (to be looked up into the `VoxelSet` primitive map)
    /// of the primitives intersected by this voxel.
    pub(crate) intersections_range: (usize, usize),
}

impl Default for Voxel {
    fn default() -> Self {
        Self {
            coords: Point::origin(),
            is_on_surface: false,
            intersections_range: (0, 0),
        }
    }
}

/// A sparse set of voxels.
///
/// It only contains voxels that are considered as "full" after a voxelization.
pub struct VoxelSet {
    /// The 3D origin of this voxel-set.
    pub origin: Point<Real>,
    /// The scale factor between the voxel integer coordinates and their
    /// actual float world-space coordinates.
    pub scale: Real,
    pub(crate) min_bb_voxels: Point<u32>,
    pub(crate) max_bb_voxels: Point<u32>,
    pub(crate) voxels: Vec<Voxel>,
    pub(crate) intersections: Arc<Vec<u32>>,
    pub(crate) primitive_classes: Arc<Vec<u32>>,
}

impl VoxelSet {
    /// Creates a new empty set of voxels.
    pub fn new() -> Self {
        Self {
            origin: Point::origin(),
            min_bb_voxels: Point::origin(),
            max_bb_voxels: Vector::repeat(1).into(),
            scale: 1.0,
            voxels: Vec::new(),
            intersections: Arc::new(Vec::new()),
            primitive_classes: Arc::new(Vec::new()),
        }
    }

    /// The volume of a single voxel of this voxel set.
    #[cfg(feature = "dim2")]
    pub fn voxel_volume(&self) -> Real {
        self.scale * self.scale
    }

    /// The volume of a single voxel of this voxel set.
    #[cfg(feature = "dim3")]
    pub fn voxel_volume(&self) -> Real {
        self.scale * self.scale * self.scale
    }

    /// Voxelizes the given shape described by its boundary:
    /// a triangle mesh (in 3D) or polyline (in 2D).
    ///
    /// # Parameters
    /// * `points` - The vertex buffer of the boundary of the shape to voxelize.
    /// * `indices` - The index buffer of the boundary of the shape to voxelize.
    /// * `resolution` - Controls the number of subdivision done along each axis. This number
    ///    is the number of subdivisions along the axis where the input shape has the largest extent.
    ///    The other dimensions will have a different automatically-determined resolution (in order to
    ///    keep the voxels cubic).
    /// * `fill_mode` - Controls what is being voxelized.
    /// * `keep_voxel_to_primitives_map` - If set to `true` a map between the voxels
    ///   and the primitives (3D triangles or 2D segments) it intersects will be computed.
    pub fn voxelize(
        points: &[Point<Real>],
        indices: &[[u32; DIM]],
        resolution: u32,
        fill_mode: FillMode,
        keep_voxel_to_primitives_map: bool,
    ) -> Self {
        VoxelizedVolume::voxelize(
            points,
            indices,
            resolution,
            fill_mode,
            keep_voxel_to_primitives_map,
        )
        .into()
    }

    /// The minimal coordinates of the integer bounding-box of the voxels in this set.
    pub fn min_bb_voxels(&self) -> Point<u32> {
        self.min_bb_voxels
    }

    /// The maximal coordinates of the integer bounding-box of the voxels in this set.
    pub fn max_bb_voxels(&self) -> Point<u32> {
        self.max_bb_voxels
    }

    /// Computes the total volume of the voxels contained by this set.
    pub fn compute_volume(&self) -> Real {
        self.voxel_volume() * self.voxels.len() as Real
    }

    fn get_voxel_point(&self, voxel: &Voxel) -> Point<Real> {
        self.get_point(na::convert(voxel.coords))
    }

    pub(crate) fn get_point(&self, voxel: Point<Real>) -> Point<Real> {
        self.origin + voxel.coords * self.scale
    }

    /// The number of voxels in this set.
    pub fn len(&self) -> usize {
        self.voxels.len()
    }

    /// The set of voxels.
    pub fn voxels(&self) -> &[Voxel] {
        &self.voxels
    }

    /// Update the bounding box of this voxel set.
    pub fn compute_bb(&mut self) {
        let num_voxels = self.voxels.len();

        if num_voxels == 0 {
            return;
        }

        self.min_bb_voxels = self.voxels[0].coords;
        self.max_bb_voxels = self.voxels[0].coords;

        for p in 0..num_voxels {
            self.min_bb_voxels = self.min_bb_voxels.inf(&self.voxels[p].coords);
            self.max_bb_voxels = self.max_bb_voxels.sup(&self.voxels[p].coords);
        }
    }

    // We have these cfg because we need to
    // use the correct return type. We could just
    // return ConvexHull but that would expose though
    // the API a type alias that isn't really worth
    // existing.
    /// Compute the convex-hull of this voxel set after cutting each voxel
    /// by the primitives (3D triangle or 2D segments) it intersects.
    ///
    /// This will panic if this `VoxelSet` was created with `keep_voxel_to_primitives_map = false`.
    #[cfg(feature = "dim2")]
    pub fn compute_exact_convex_hull(
        &self,
        points: &[Point<Real>],
        indices: &[[u32; DIM]],
    ) -> Vec<Point<Real>> {
        self.do_compute_exact_convex_hull(points, indices)
    }

    /// Compute the convex-hull of this voxel set after cutting each voxel
    /// by the primitives (3D triangle or 2D segments) it intersects.
    ///
    /// This will panic if this `VoxelSet` was created with `keep_voxel_to_primitives_map = false`.
    #[cfg(feature = "dim3")]
    pub fn compute_exact_convex_hull(
        &self,
        points: &[Point<Real>],
        indices: &[[u32; DIM]],
    ) -> (Vec<Point<Real>>, Vec<[u32; DIM]>) {
        self.do_compute_exact_convex_hull(points, indices)
    }

    fn do_compute_exact_convex_hull(
        &self,
        points: &[Point<Real>],
        indices: &[[u32; DIM]],
    ) -> ConvexHull {
        assert!(!self.intersections.is_empty(),
                "Cannot compute exact convex hull without voxel-to-primitives-map. Consider passing voxel_to_primitives_map = true to the voxelizer.");
        let mut surface_points = Vec::new();
        #[cfg(feature = "dim3")]
        let (mut polygon, mut workspace) = (Vec::new(), Vec::new());
        let mut pushed_points = vec![false; points.len()];

        // Grab all the points.
        for voxel in self.voxels.iter().filter(|v| v.is_on_surface) {
            let intersections =
                &self.intersections[voxel.intersections_range.0..voxel.intersections_range.1];
            for prim_id in intersections {
                let ia = indices[*prim_id as usize][0] as usize;
                let ib = indices[*prim_id as usize][1] as usize;
                #[cfg(feature = "dim3")]
                let ic = indices[*prim_id as usize][2] as usize;

                // If the primitives have been classified by VHACD, we know that:
                // - A class equal to Some(u32::MAX) means that the primitives intersects multiple
                //   convex parts, so we need to split it.
                // - A class equal to None means that we did not compute any classes (so we
                //   must assume that each triangle have to be split since it may intersect
                //   multiple parts.
                // - A class different from `None` and `Some(u32::MAX)` means that the triangle is
                //   included in only one convex part. So instead of cutting it, just push the whole
                //   triangle once.
                let prim_class = self.primitive_classes.get(*prim_id as usize).copied();
                if prim_class == Some(u32::MAX) || prim_class == None {
                    let aabb_center =
                        self.origin + voxel.coords.coords.map(|k| k as Real) * self.scale;
                    let aabb =
                        AABB::from_half_extents(aabb_center, Vector::repeat(self.scale / 2.0));

                    #[cfg(feature = "dim2")]
                    if let Some(seg) = aabb.clip_segment(&points[ia], &points[ib]) {
                        surface_points.push(seg.a);
                        surface_points.push(seg.b);
                    }

                    #[cfg(feature = "dim3")]
                    {
                        polygon.clear();
                        polygon.extend_from_slice(&[points[ia], points[ib], points[ic]]);
                        aabb.clip_polygon_with_workspace(&mut polygon, &mut workspace);
                        surface_points.append(&mut polygon);
                    }
                } else {
                    // We know this triangle is only contained by
                    // one voxel set, i.e., `self`. So we don't
                    // need to cut it.
                    //
                    // Because one triangle may intersect multiple voxels contained by
                    // the same convex part, we only push vertices we have not pushed
                    // so far in order to avoid some useless duplicate points (duplicate
                    // points are OK as far as convex hull computation is concerned, but
                    // they imply some redundant computations).
                    let mut push_pt = |i: usize| {
                        if !pushed_points[i] {
                            surface_points.push(points[i]);
                            pushed_points[i] = true;
                        }
                    };

                    push_pt(ia);
                    push_pt(ib);
                    #[cfg(feature = "dim3")]
                    push_pt(ic);
                }
            }

            if intersections.is_empty() {
                self.map_voxel_points(voxel, |p| surface_points.push(p));
            }
        }

        // Compute the convex-hull.
        convex_hull(&surface_points)
    }

    /// Computes the intersections between all the voxels of this voxel set,
    /// and all the primitives (triangle or segments) it intersected (as per
    /// the voxel-to-primitives-map computed during voxelization).
    ///
    /// Panics if the voxelization was performed without setting the parameter
    /// `voxel_to_primitives_map = true`.
    pub fn compute_primitive_intersections(
        &self,
        points: &[Point<Real>],
        indices: &[[u32; DIM]],
    ) -> Vec<Point<Real>> {
        assert!(!self.intersections.is_empty(),
                "Cannot compute primitive intersections voxel-to-primitives-map. Consider passing voxel_to_primitives_map = true to the voxelizer.");
        let mut surface_points = Vec::new();
        #[cfg(feature = "dim3")]
        let (mut polygon, mut workspace) = (Vec::new(), Vec::new());

        // Grab all the points.
        for voxel in self.voxels.iter().filter(|v| v.is_on_surface) {
            let intersections =
                &self.intersections[voxel.intersections_range.0..voxel.intersections_range.1];
            for prim_id in intersections {
                let aabb_center = self.origin + voxel.coords.coords.map(|k| k as Real) * self.scale;
                let aabb = AABB::from_half_extents(aabb_center, Vector::repeat(self.scale / 2.0));

                let pa = points[indices[*prim_id as usize][0] as usize];
                let pb = points[indices[*prim_id as usize][1] as usize];
                #[cfg(feature = "dim3")]
                let pc = points[indices[*prim_id as usize][2] as usize];

                #[cfg(feature = "dim2")]
                if let Some(seg) = aabb.clip_segment(&pa, &pb) {
                    surface_points.push(seg.a);
                    surface_points.push(seg.b);
                }

                #[cfg(feature = "dim3")]
                {
                    workspace.clear();
                    polygon.clear();
                    polygon.extend_from_slice(&[pa, pb, pc]);
                    aabb.clip_polygon_with_workspace(&mut polygon, &mut workspace);

                    if polygon.len() > 2 {
                        for i in 1..polygon.len() - 1 {
                            surface_points.push(polygon[0]);
                            surface_points.push(polygon[i]);
                            surface_points.push(polygon[i + 1]);
                        }
                    }
                }
            }
        }

        surface_points
    }

    /// Compute the convex-hull of the voxels in this set.
    ///
    /// # Parameters
    /// * `sampling` - The convex-hull computation will ignore `sampling` voxels at
    ///   regular intervals. Useful to save some computation times if an exact result isn't need.
    ///   Use `0` to make sure no voxel is being ignored.
    #[cfg(feature = "dim2")]
    pub fn compute_convex_hull(&self, sampling: u32) -> Vec<Point<Real>> {
        let mut points = Vec::new();

        // Grab all the points.
        for voxel in self
            .voxels
            .iter()
            .filter(|v| v.is_on_surface)
            .step_by(sampling as usize)
        {
            self.map_voxel_points(voxel, |p| points.push(p));
        }

        // Compute the convex-hull.
        convex_hull(&points)
    }

    /// Compute the convex-hull of the voxels in this set.
    ///
    /// # Parameters
    /// * `sampling` - The convex-hull computation will ignore `sampling` voxels at
    ///   regular intervals. Useful to save some computation times if an exact result isn't need.
    ///   Use `0` to make sure no voxel is being ignored.
    #[cfg(feature = "dim3")]
    pub fn compute_convex_hull(&self, sampling: u32) -> (Vec<Point<Real>>, Vec<[u32; DIM]>) {
        let mut points = Vec::new();

        // Grab all the points.
        for voxel in self
            .voxels
            .iter()
            .filter(|v| v.is_on_surface)
            .step_by(sampling as usize)
        {
            self.map_voxel_points(voxel, |p| points.push(p));
        }

        // Compute the convex-hull.
        convex_hull(&points)
    }

    /// Gets the vertices of the given voxel.
    fn map_voxel_points(&self, voxel: &Voxel, mut f: impl FnMut(Point<Real>)) {
        let ijk = voxel.coords.coords.map(|e| e as Real);

        #[cfg(feature = "dim2")]
        let shifts = [
            Vector::new(-0.5, -0.5),
            Vector::new(0.5, -0.5),
            Vector::new(0.5, 0.5),
            Vector::new(-0.5, 0.5),
        ];

        #[cfg(feature = "dim3")]
        let shifts = [
            Vector::new(-0.5, -0.5, -0.5),
            Vector::new(0.5, -0.5, -0.5),
            Vector::new(0.5, 0.5, -0.5),
            Vector::new(-0.5, 0.5, -0.5),
            Vector::new(-0.5, -0.5, 0.5),
            Vector::new(0.5, -0.5, 0.5),
            Vector::new(0.5, 0.5, 0.5),
            Vector::new(-0.5, 0.5, 0.5),
        ];

        for shift in &shifts {
            f(self.origin + (ijk + *shift) * self.scale)
        }
    }

    pub(crate) fn intersect(
        &self,
        plane: &CutPlane,
        positive_pts: &mut Vec<Point<Real>>,
        negative_pts: &mut Vec<Point<Real>>,
        sampling: u32,
    ) {
        let num_voxels = self.voxels.len();

        if num_voxels == 0 {
            return;
        }

        let d0 = self.scale;
        let mut sp = 0;
        let mut sn = 0;

        for v in 0..num_voxels {
            let voxel = self.voxels[v];
            let pt = self.get_voxel_point(&voxel);
            let d = plane.abc.dot(&pt.coords) + plane.d;

            // if      (d >= 0.0 && d <= d0) positive_pts.push(pt);
            // else if (d < 0.0 && -d <= d0) negative_pts.push(pt);

            if d >= 0.0 {
                if d <= d0 {
                    self.map_voxel_points(&voxel, |p| positive_pts.push(p));
                } else {
                    sp += 1;

                    if sp == sampling {
                        self.map_voxel_points(&voxel, |p| positive_pts.push(p));
                        sp = 0;
                    }
                }
            } else {
                if -d <= d0 {
                    self.map_voxel_points(&voxel, |p| negative_pts.push(p));
                } else {
                    sn += 1;
                    if sn == sampling {
                        self.map_voxel_points(&voxel, |p| negative_pts.push(p));
                        sn = 0;
                    }
                }
            }
        }
    }

    // Returns (negative_volume, positive_volume)
    pub(crate) fn compute_clipped_volumes(&self, plane: &CutPlane) -> (Real, Real) {
        if self.voxels.is_empty() {
            return (0.0, 0.0);
        }

        let mut num_positive_voxels = 0;

        for voxel in &self.voxels {
            let pt = self.get_voxel_point(voxel);
            let d = plane.abc.dot(&pt.coords) + plane.d;
            num_positive_voxels += (d >= 0.0) as usize;
        }

        let num_negative_voxels = self.voxels.len() - num_positive_voxels;
        let positive_volume = self.voxel_volume() * (num_positive_voxels as Real);
        let negative_volume = self.voxel_volume() * (num_negative_voxels as Real);

        (negative_volume, positive_volume)
    }

    // Set `on_surf` such that it contains only the voxel on surface contained by `self`.
    pub(crate) fn select_on_surface(&self, on_surf: &mut VoxelSet) {
        on_surf.origin = self.origin;
        on_surf.voxels.clear();
        on_surf.scale = self.scale;

        for voxel in &self.voxels {
            if voxel.is_on_surface {
                on_surf.voxels.push(*voxel);
            }
        }
    }

    /// Splits this voxel set into two parts, depending on where the voxel center lies wrt. the given plane.
    pub(crate) fn clip(
        &self,
        plane: &CutPlane,
        positive_part: &mut VoxelSet,
        negative_part: &mut VoxelSet,
    ) {
        let num_voxels = self.voxels.len();

        if num_voxels == 0 {
            return;
        }

        negative_part.origin = self.origin;
        negative_part.voxels.clear();
        negative_part.voxels.reserve(num_voxels);
        negative_part.scale = self.scale;

        positive_part.origin = self.origin;
        positive_part.voxels.clear();
        positive_part.voxels.reserve(num_voxels);
        positive_part.scale = self.scale;

        let d0 = self.scale;

        for v in 0..num_voxels {
            let mut voxel = self.voxels[v];
            let pt = self.get_voxel_point(&voxel);
            let d = plane.abc.dot(&pt.coords) + plane.d;

            if d >= 0.0 {
                if voxel.is_on_surface || d <= d0 {
                    voxel.is_on_surface = true;
                    positive_part.voxels.push(voxel);
                } else {
                    positive_part.voxels.push(voxel);
                }
            } else {
                if voxel.is_on_surface || -d <= d0 {
                    voxel.is_on_surface = true;
                    negative_part.voxels.push(voxel);
                } else {
                    negative_part.voxels.push(voxel);
                }
            }
        }
    }

    /// Convert `self` into a mesh, including only the voxels on the surface or only the voxel
    /// inside of the volume.
    #[cfg(feature = "dim3")]
    pub fn to_trimesh(
        &self,
        base_index: u32,
        is_on_surface: bool,
    ) -> (Vec<Point<Real>>, Vec<[u32; DIM]>) {
        let mut vertices = Vec::new();
        let mut indices = Vec::new();

        for voxel in &self.voxels {
            if voxel.is_on_surface == is_on_surface {
                self.map_voxel_points(voxel, |p| vertices.push(p));

                indices.push([base_index + 0, base_index + 2, base_index + 1]);
                indices.push([base_index + 0, base_index + 3, base_index + 2]);
                indices.push([base_index + 4, base_index + 5, base_index + 6]);
                indices.push([base_index + 4, base_index + 6, base_index + 7]);
                indices.push([base_index + 7, base_index + 6, base_index + 2]);
                indices.push([base_index + 7, base_index + 2, base_index + 3]);
                indices.push([base_index + 4, base_index + 1, base_index + 5]);
                indices.push([base_index + 4, base_index + 0, base_index + 1]);
                indices.push([base_index + 6, base_index + 5, base_index + 1]);
                indices.push([base_index + 6, base_index + 1, base_index + 2]);
                indices.push([base_index + 7, base_index + 0, base_index + 4]);
                indices.push([base_index + 7, base_index + 3, base_index + 0]);
            }
        }

        (vertices, indices)
    }

    pub(crate) fn compute_principal_axes(&self) -> Vector<Real> {
        let num_voxels = self.voxels.len();
        if num_voxels == 0 {
            return Vector::zeros();
        }

        // TODO: find a way to reuse crate::utils::cov?
        // The difficulty being that we need to iterate through the set of
        // points twice. So passing an iterator to crate::utils::cov
        // isn't really possible.
        let mut center = Point::origin();
        let denom = 1.0 / (num_voxels as Real);

        for voxel in &self.voxels {
            center += voxel.coords.map(|e| e as Real).coords * denom;
        }

        let mut cov_mat = Matrix::zeros();
        for voxel in &self.voxels {
            let xyz = voxel.coords.map(|e| e as Real) - center;
            cov_mat.syger(denom, &xyz, &xyz, 1.0);
        }

        cov_mat.symmetric_eigenvalues()
    }
}

#[cfg(feature = "dim2")]
fn convex_hull(vertices: &[Point<Real>]) -> Vec<Point<Real>> {
    if vertices.len() > 1 {
        crate::transformation::convex_hull(vertices)
    } else {
        Vec::new()
    }
}

#[cfg(feature = "dim3")]
fn convex_hull(vertices: &[Point<Real>]) -> (Vec<Point<Real>>, Vec<[u32; DIM]>) {
    if vertices.len() > 2 {
        crate::transformation::convex_hull(vertices)
    } else {
        (Vec::new(), Vec::new())
    }
}
