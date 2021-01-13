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

use super::{FillMode, VoxelValue, VoxelizedVolume};
use crate::bounding_volume::AABB;
use crate::math::Real;
use crate::na::{Isometry3, SymmetricEigen};
use crate::query;
use crate::shape::Triangle;
use crate::transformation::vhacd::Plane;
use na::{Matrix3, Point3, Vector3};

#[derive(Copy, Clone, Debug)]
pub struct Voxel {
    pub coords: Point3<u32>,
    pub data: VoxelValue,
}

impl Default for Voxel {
    fn default() -> Self {
        Self {
            coords: Point3::origin(),
            data: VoxelValue::PrimitiveUndefined,
        }
    }
}

pub struct VoxelSet {
    pub min_bb: Point3<Real>,
    pub(crate) min_bb_voxels: Point3<u32>,
    pub(crate) max_bb_voxels: Point3<u32>,
    pub(crate) min_bb_pts: Point3<Real>,
    pub(crate) max_bb_pts: Point3<Real>,
    pub(crate) barycenter: Point3<u32>,
    pub(crate) barycenter_pca: Point3<Real>,
    pub scale: Real,
    pub(crate) unit_volume: Real,
    pub(crate) num_voxels_on_surface: u32,
    pub(crate) num_voxels_inside_surface: u32,
    pub(crate) voxels: Vec<Voxel>,
    pub(crate) eigenvalues: Vector3<Real>,
}

impl VoxelSet {
    pub fn new() -> Self {
        Self {
            min_bb: Point3::origin(),
            min_bb_voxels: Point3::origin(),
            max_bb_voxels: Point3::new(1, 1, 1),
            min_bb_pts: Point3::origin(),
            max_bb_pts: Point3::new(1.0, 1.0, 1.0),
            barycenter: Point3::origin(),
            barycenter_pca: Point3::origin(),
            scale: 1.0,
            unit_volume: 1.0,
            num_voxels_on_surface: 0,
            num_voxels_inside_surface: 0,
            voxels: Vec::new(),
            eigenvalues: Vector3::zeros(),
        }
    }

    pub fn voxelize(
        transform: &Isometry3<Real>,
        points: &[Point3<Real>],
        triangles: &[Point3<u32>],
        dim: u32,
        fill_mode: FillMode,
    ) -> Self {
        VoxelizedVolume::voxelize(transform, points, triangles, dim, fill_mode).into()
    }

    pub fn min_bb_voxels(&self) -> Point3<u32> {
        self.min_bb_voxels
    }

    pub fn max_bb_voxels(&self) -> Point3<u32> {
        self.max_bb_voxels
    }

    pub fn eigenvalues(&self) -> Vector3<Real> {
        self.eigenvalues
    }

    pub fn compute_volume(&self) -> Real {
        self.unit_volume * self.voxels.len() as Real
    }

    fn get_voxel_point(&self, voxel: &Voxel) -> Point3<Real> {
        self.get_point(na::convert(voxel.coords))
    }

    pub fn get_point(&self, voxel: Point3<Real>) -> Point3<Real> {
        self.min_bb + voxel.coords * self.scale
    }

    pub fn len(&self) -> usize {
        self.voxels.len()
    }

    /// Update the bounding box of this voxel set.
    pub fn compute_bb(&mut self) {
        let num_voxels = self.voxels.len();

        if num_voxels == 0 {
            return;
        }

        self.min_bb_voxels = self.voxels[0].coords;
        self.max_bb_voxels = self.voxels[0].coords;
        let mut bary = self.voxels[0].coords;

        for p in 0..num_voxels {
            bary += self.voxels[p].coords.coords;
            self.min_bb_voxels = self.min_bb_voxels.inf(&self.voxels[p].coords);
            self.max_bb_voxels = self.max_bb_voxels.sup(&self.voxels[p].coords);
        }

        let bary = bary.coords.map(|e| e as Real / num_voxels as Real);
        self.min_bb_pts = self.min_bb + self.min_bb_voxels.coords.map(|e| e as Real * self.scale);
        self.max_bb_pts = self.min_bb + self.max_bb_voxels.coords.map(|e| e as Real * self.scale);
    }

    pub fn compute_convex_hull(&self, sampling: u32) -> (Vec<Point3<Real>>, Vec<Point3<u32>>) {
        let mut points = Vec::new();

        // Grab all the points.
        for voxel in self
            .voxels
            .iter()
            .filter(|v| v.data == VoxelValue::PrimitiveOnSurface)
            .step_by(sampling as usize)
        {
            self.map_voxel_points(voxel, |p| points.push(p));
        }

        // Compute the convex-hull.
        convex_hull(&points)
    }

    /// Gets the vertices of the given voxel.
    fn map_voxel_points(&self, voxel: &Voxel, mut f: impl FnMut(Point3<Real>)) {
        let ijk = voxel.coords.coords.map(|e| e as Real);

        let shifts = [
            Vector3::new(-0.5, -0.5, -0.5),
            Vector3::new(0.5, -0.5, -0.5),
            Vector3::new(0.5, 0.5, -0.5),
            Vector3::new(-0.5, 0.5, -0.5),
            Vector3::new(-0.5, -0.5, 0.5),
            Vector3::new(0.5, -0.5, 0.5),
            Vector3::new(0.5, 0.5, 0.5),
            Vector3::new(-0.5, 0.5, 0.5),
        ];

        for l in 0..8 {
            f(self.min_bb + (ijk + shifts[l]) * self.scale)
        }
    }

    pub fn intersect(
        &self,
        plane: &Plane,
        positive_pts: &mut Vec<Point3<Real>>,
        negative_pts: &mut Vec<Point3<Real>>,
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
    pub fn compute_clipped_volumes(&self, plane: &Plane) -> (Real, Real) {
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
        let positive_volume = self.unit_volume * (num_positive_voxels as Real);
        let negative_volume = self.unit_volume * (num_negative_voxels as Real);

        (negative_volume, positive_volume)
    }

    // Set `on_surf` such that it contains only the voxel on surface contained by `self`.
    pub fn select_on_surface(&self, on_surf: &mut VoxelSet) {
        on_surf.min_bb = self.min_bb;
        on_surf.voxels.clear();
        on_surf.scale = self.scale;
        on_surf.unit_volume = self.unit_volume;
        on_surf.num_voxels_on_surface = 0;
        on_surf.num_voxels_inside_surface = 0;

        for voxel in &self.voxels {
            if voxel.data == VoxelValue::PrimitiveOnSurface {
                on_surf.voxels.push(*voxel);
                on_surf.num_voxels_on_surface += 1;
            }
        }
    }

    /// Splits this voxel set into two parts, depending on where the voxel center lies wrt. the given plane.
    pub fn clip(&self, plane: &Plane, positive_part: &mut VoxelSet, negative_part: &mut VoxelSet) {
        let num_voxels = self.voxels.len();

        if num_voxels == 0 {
            return;
        }

        negative_part.min_bb = self.min_bb;
        negative_part.voxels.clear();
        negative_part.voxels.reserve(num_voxels);
        negative_part.scale = self.scale;
        negative_part.unit_volume = self.unit_volume;
        negative_part.num_voxels_on_surface = 0;
        negative_part.num_voxels_inside_surface = 0;

        positive_part.min_bb = self.min_bb;
        positive_part.voxels.clear();
        positive_part.voxels.reserve(num_voxels);
        positive_part.scale = self.scale;
        positive_part.unit_volume = self.unit_volume;
        positive_part.num_voxels_on_surface = 0;
        positive_part.num_voxels_inside_surface = 0;

        let d0 = self.scale;

        for v in 0..num_voxels {
            let mut voxel = self.voxels[v];
            let pt = self.get_voxel_point(&voxel);
            let d = plane.abc.dot(&pt.coords) + plane.d;

            if d >= 0.0 {
                if voxel.data == VoxelValue::PrimitiveOnSurface || d <= d0 {
                    voxel.data = VoxelValue::PrimitiveOnSurface;
                    positive_part.voxels.push(voxel);
                    positive_part.num_voxels_on_surface += 1;
                } else {
                    positive_part.voxels.push(voxel);
                    positive_part.num_voxels_inside_surface += 1;
                }
            } else {
                if voxel.data == VoxelValue::PrimitiveOnSurface || -d <= d0 {
                    voxel.data = VoxelValue::PrimitiveOnSurface;
                    negative_part.voxels.push(voxel);
                    negative_part.num_voxels_on_surface += 1;
                } else {
                    negative_part.voxels.push(voxel);
                    negative_part.num_voxels_inside_surface += 1;
                }
            }
        }
    }

    /// Convert this voxelset into a mesh, including only the voxels with the given value.
    fn to_trimesh(
        &self,
        base_index: u32,
        value: VoxelValue,
    ) -> (Vec<Point3<Real>>, Vec<Point3<u32>>) {
        let mut vertices = Vec::new();
        let mut indices = Vec::new();

        for voxel in &self.voxels {
            if voxel.data == value {
                self.map_voxel_points(voxel, |p| vertices.push(p));

                indices.push(Point3::new(base_index + 0, base_index + 2, base_index + 1));
                indices.push(Point3::new(base_index + 0, base_index + 3, base_index + 2));
                indices.push(Point3::new(base_index + 4, base_index + 5, base_index + 6));
                indices.push(Point3::new(base_index + 4, base_index + 6, base_index + 7));
                indices.push(Point3::new(base_index + 7, base_index + 6, base_index + 2));
                indices.push(Point3::new(base_index + 7, base_index + 2, base_index + 3));
                indices.push(Point3::new(base_index + 4, base_index + 1, base_index + 5));
                indices.push(Point3::new(base_index + 4, base_index + 0, base_index + 1));
                indices.push(Point3::new(base_index + 6, base_index + 5, base_index + 1));
                indices.push(Point3::new(base_index + 6, base_index + 1, base_index + 2));
                indices.push(Point3::new(base_index + 7, base_index + 0, base_index + 4));
                indices.push(Point3::new(base_index + 7, base_index + 3, base_index + 0));
            }
        }

        (vertices, indices)
    }

    pub fn compute_principal_axes(&mut self) {
        let num_voxels = self.voxels.len();
        if num_voxels == 0 {
            return;
        }

        self.barycenter_pca = Point3::origin();
        let denom = 1.0 / (num_voxels as Real);

        for voxel in &self.voxels {
            self.barycenter_pca += voxel.coords.map(|e| e as Real).coords * denom;
        }

        let mut cov_mat = Matrix3::zeros();

        for voxel in &self.voxels {
            let xyz = voxel.coords.map(|e| e as Real) - self.barycenter_pca;
            cov_mat.syger(denom, &xyz, &xyz, 1.0);
        }

        self.eigenvalues = cov_mat.symmetric_eigenvalues();
    }
}

fn convex_hull(vertices: &[Point3<Real>]) -> (Vec<Point3<Real>>, Vec<Point3<u32>>) {
    if vertices.len() > 2 {
        crate::transformation::convex_hull(vertices)
    } else {
        (Vec::new(), Vec::new())
    }
}
