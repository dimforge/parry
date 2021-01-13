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

use crate::bounding_volume::AABB;
use crate::math::Real;
use crate::na::Isometry3;
use crate::query;
use crate::shape::Triangle;
use crate::transformation::voxelization::{Voxel, VoxelSet};
use na::{Point3, Vector3};

#[derive(Copy, Clone, PartialEq, Eq)]
pub enum FillMode {
    SurfaceOnly,
    FloodFill,
    // RaycastFill
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum VoxelValue {
    PrimitiveUndefined,
    PrimitiveOutsideSurfaceToWalk,
    PrimitiveOutsideSurface,
    PrimitiveInsideSurface,
    PrimitiveOnSurface,
}

pub struct VoxelizedVolume {
    origin: Point3<Real>,
    scale: Real,
    resolution: Point3<u32>,
    num_voxels_on_surface: u32,
    num_voxels_inside_surface: u32,
    num_voxels_outside_surface: u32,
    data: Vec<VoxelValue>,
}

impl VoxelizedVolume {
    pub fn voxelize(
        transform: &Isometry3<Real>,
        points: &[Point3<Real>],
        triangles: &[Point3<u32>],
        resolution: u32,
        fill_mode: FillMode,
    ) -> Self {
        let mut result = VoxelizedVolume {
            resolution: Point3::origin(),
            origin: Point3::origin(),
            num_voxels_on_surface: 0,
            num_voxels_inside_surface: 0,
            num_voxels_outside_surface: 0,
            scale: 1.0,
            data: Vec::new(),
        };

        if points.is_empty() {
            return result;
        }

        let aabb = crate::bounding_volume::point_cloud_aabb(transform, points);
        result.origin = aabb.mins;

        let d = aabb.maxs - aabb.mins;
        let r;

        if d[0] > d[1] && d[0] > d[2] {
            r = d[0];
            result.resolution[0] = resolution;
            result.resolution[1] = 2 + (resolution as Real * d[1] / d[0]) as u32;
            result.resolution[2] = 2 + (resolution as Real * d[2] / d[0]) as u32;
        } else if d[1] > d[0] && d[1] > d[2] {
            r = d[1];
            result.resolution[1] = resolution;
            result.resolution[0] = 2 + (resolution as Real * d[0] / d[1]) as u32;
            result.resolution[2] = 2 + (resolution as Real * d[2] / d[1]) as u32;
        } else {
            r = d[2];
            result.resolution[2] = resolution;
            result.resolution[0] = 2 + (resolution as Real * d[0] / d[2]) as u32;
            result.resolution[1] = 2 + (resolution as Real * d[1] / d[2]) as u32;
        }

        result.scale = r / (resolution as Real - 1.0);
        let inv_scale = (resolution as Real - 1.0) / r;
        result.allocate();

        result.num_voxels_on_surface = 0;
        result.num_voxels_on_surface = 0;
        result.num_voxels_inside_surface = 0;
        result.num_voxels_outside_surface = 0;

        let mut tri_pts = [Point3::origin(); 3];
        let box_half_size = Vector3::repeat(0.5);
        let mut ijk0 = Vector3::new(0u32, 0, 0);
        let mut ijk1 = Vector3::new(0u32, 0, 0);

        for tri in triangles {
            // Find the range of voxels potentially intersecting the triangle.
            for c in 0..3 {
                let pt = transform * points[tri[c] as usize];
                tri_pts[c] = (pt - result.origin.coords) * inv_scale;

                let i = (tri_pts[c].x + 0.5) as u32;
                let j = (tri_pts[c].y + 0.5) as u32;
                let k = (tri_pts[c].z + 0.5) as u32;

                assert!(
                    i < result.resolution[0]
                        && j < result.resolution[1]
                        && k < result.resolution[2]
                );

                let ijk = Vector3::new(i, j, k);

                if c == 0 {
                    ijk0 = ijk;
                    ijk1 = ijk;
                } else {
                    ijk0 = ijk0.inf(&ijk);
                    ijk1 = ijk1.sup(&ijk);
                }
            }

            ijk0.apply(|e| e.saturating_sub(1));
            ijk1 = ijk1.map(|e| e + 1).inf(&result.resolution.coords);

            // Determine exactly what voxel intersect the triangle.
            for i in ijk0.x..ijk1.x {
                for j in ijk0.y..ijk1.y {
                    for k in ijk0.z..ijk1.z {
                        let value = result.voxel_mut(i, j, k);

                        if *value == VoxelValue::PrimitiveUndefined {
                            let triangle = Triangle::from(tri_pts);
                            let aabb = AABB::from_half_extents(
                                Point3::new(i as Real, j as Real, k as Real),
                                box_half_size,
                            );

                            if query::details::intersection_test_aabb_triangle(&aabb, &triangle) {
                                *value = VoxelValue::PrimitiveOnSurface;
                                result.num_voxels_on_surface += 1;
                            }
                        }
                    }
                }
            }
        }

        match fill_mode {
            FillMode::SurfaceOnly => {
                for value in &mut result.data {
                    if *value != VoxelValue::PrimitiveOnSurface {
                        *value = VoxelValue::PrimitiveOutsideSurface
                    }
                }
            }
            FillMode::FloodFill => {
                result.mark_outside_surface(0, 0, 0, result.resolution[0], result.resolution[1], 1);
                result.mark_outside_surface(
                    0,
                    0,
                    result.resolution[2] - 1,
                    result.resolution[0],
                    result.resolution[1],
                    result.resolution[2],
                );
                result.mark_outside_surface(0, 0, 0, result.resolution[0], 1, result.resolution[2]);
                result.mark_outside_surface(
                    0,
                    result.resolution[1] - 1,
                    0,
                    result.resolution[0],
                    result.resolution[1],
                    result.resolution[2],
                );
                result.mark_outside_surface(0, 0, 0, 1, result.resolution[1], result.resolution[2]);
                result.mark_outside_surface(
                    result.resolution[0] - 1,
                    0,
                    0,
                    result.resolution[0],
                    result.resolution[1],
                    result.resolution[2],
                );
                result.fill_outside_surface();
                result.fill_inside_surface();
            }
        }

        result
    }

    pub fn resolution(&self) -> Point3<u32> {
        self.resolution
    }

    pub fn scale(&self) -> Real {
        self.scale
    }

    fn allocate(&mut self) {
        let len = self.resolution[0] * self.resolution[1] * self.resolution[2];
        self.data
            .resize(len as usize, VoxelValue::PrimitiveUndefined);
    }

    fn voxel_index(&self, i: u32, j: u32, k: u32) -> u32 {
        i + j * self.resolution[0] + k * self.resolution[0] * self.resolution[1]
    }

    fn voxel_mut(&mut self, i: u32, j: u32, k: u32) -> &mut VoxelValue {
        let idx = self.voxel_index(i, j, k);
        &mut self.data[idx as usize]
    }

    pub fn voxel(&self, i: u32, j: u32, k: u32) -> VoxelValue {
        let idx = self.voxel_index(i, j, k);
        self.data[idx as usize]
    }

    pub fn num_voxels_on_surface(&self) -> u32 {
        self.num_voxels_on_surface
    }

    pub fn num_voxels_inside_surface(&self) -> u32 {
        self.num_voxels_inside_surface
    }

    pub fn num_voxels_outside_surface(&self) -> u32 {
        self.num_voxels_outside_surface
    }
    /// Mark all the PrimitiveUndefined voxels within the given bounds as PrimitiveOutsideSurfaceToWalk.
    fn mark_outside_surface(&mut self, i0: u32, j0: u32, k0: u32, i1: u32, j1: u32, k1: u32) {
        for i in i0..i1 {
            for j in j0..j1 {
                for k in k0..k1 {
                    let v = self.voxel_mut(i, j, k);

                    if *v == VoxelValue::PrimitiveUndefined {
                        *v = VoxelValue::PrimitiveOutsideSurfaceToWalk;
                    }
                }
            }
        }
    }

    fn walk_forward(
        start: isize,
        end: isize,
        mut ptr: isize,
        out: &mut [VoxelValue],
        stride: isize,
        max_distance: isize,
    ) {
        let mut i = start;
        let mut count = 0;

        while count < max_distance && i < end && out[ptr as usize] == VoxelValue::PrimitiveUndefined
        {
            out[ptr as usize] = VoxelValue::PrimitiveOutsideSurfaceToWalk;
            i += 1;
            ptr += stride;
            count += 1;
        }
    }

    fn walk_backward(
        start: isize,
        end: isize,
        mut ptr: isize,
        out: &mut [VoxelValue],
        stride: isize,
        max_distance: isize,
    ) {
        let mut i = start;
        let mut count = 0;

        while count < max_distance
            && i >= end
            && out[ptr as usize] == VoxelValue::PrimitiveUndefined
        {
            out[ptr as usize] = VoxelValue::PrimitiveOutsideSurfaceToWalk;
            i -= 1;
            ptr -= stride;
            count += 1;
        }
    }

    fn fill_outside_surface(&mut self) {
        let mut voxels_walked;
        let i0 = self.resolution[0];
        let j0 = self.resolution[1];
        let k0 = self.resolution[2];

        // Avoid striding too far in each direction to stay in L1 cache as much as possible.
        // The cache size required for the walk is roughly (4 * walk_distance * 64) since
        // the k direction doesn't count as it's walking byte per byte directly in a cache lines.
        // ~16k is required for a walk distance of 64 in each directions.
        let walk_distance = 64;

        // using the stride directly instead of calling get_voxel for each iterations saves
        // a lot of multiplications and pipeline stalls due to data dependencies on imul.
        let istride = self.voxel_index(1, 0, 0) as isize - self.voxel_index(0, 0, 0) as isize;
        let jstride = self.voxel_index(0, 1, 0) as isize - self.voxel_index(0, 0, 0) as isize;
        let kstride = self.voxel_index(0, 0, 1) as isize - self.voxel_index(0, 0, 0) as isize;

        // It might seem counter intuitive to go over the whole voxel range multiple times
        // but since we do the run in memory order, it leaves us with far fewer cache misses
        // than a BFS algorithm and it has the additional benefit of not requiring us to
        // store and manipulate a fifo for recursion that might become huge when the number
        // of voxels is large.
        // This will outperform the BFS algorithm by several orders of magnitude in practice.
        loop {
            voxels_walked = 0;

            for i in 0..i0 {
                for j in 0..j0 {
                    for k in 0..k0 {
                        let idx = self.voxel_index(i, j, k) as isize;
                        let voxel = self.voxel_mut(i, j, k);

                        if *voxel == VoxelValue::PrimitiveOutsideSurfaceToWalk {
                            voxels_walked += 1;
                            *voxel = VoxelValue::PrimitiveOutsideSurface;

                            // walk in each direction to mark other voxel that should be walked.
                            // this will generate a 3d pattern that will help the overall
                            // algorithm converge faster while remaining cache friendly.
                            Self::walk_forward(
                                k as isize + 1,
                                k0 as isize,
                                idx + kstride,
                                &mut self.data,
                                kstride,
                                walk_distance,
                            );
                            Self::walk_backward(
                                k as isize - 1,
                                0,
                                idx - kstride,
                                &mut self.data,
                                kstride,
                                walk_distance,
                            );

                            Self::walk_forward(
                                j as isize + 1,
                                j0 as isize,
                                idx + jstride,
                                &mut self.data,
                                jstride,
                                walk_distance,
                            );
                            Self::walk_backward(
                                j as isize - 1,
                                0,
                                idx - jstride,
                                &mut self.data,
                                jstride,
                                walk_distance,
                            );

                            Self::walk_forward(
                                (i + 1) as isize,
                                i0 as isize,
                                idx + istride,
                                &mut self.data,
                                istride,
                                walk_distance,
                            );
                            Self::walk_backward(
                                i as isize - 1,
                                0,
                                idx - istride,
                                &mut self.data,
                                istride,
                                walk_distance,
                            );
                        }
                    }
                }
            }

            self.num_voxels_outside_surface += voxels_walked;

            if voxels_walked == 0 {
                break;
            }
        }
    }

    fn fill_inside_surface(&mut self) {
        for i in 0..self.resolution.x {
            for j in 0..self.resolution.y {
                for k in 0..self.resolution.z {
                    let v = self.voxel_mut(i, j, k);
                    if *v == VoxelValue::PrimitiveUndefined {
                        *v = VoxelValue::PrimitiveInsideSurface;
                        self.num_voxels_inside_surface += 1;
                    }
                }
            }
        }
    }

    pub fn to_trimesh(&self, value: VoxelValue) -> (Vec<Point3<Real>>, Vec<Point3<u32>>) {
        let mut vertices = Vec::new();
        let mut indices = Vec::new();

        for i in 0..self.resolution.x {
            for j in 0..self.resolution.y {
                for k in 0..self.resolution.z {
                    let voxel = self.voxel(i, j, k);

                    if voxel == value {
                        let ijk = Vector3::new(i as Real, j as Real, k as Real);

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

                        for shift in &shifts {
                            vertices.push(self.origin + (ijk + shift) * self.scale);
                        }

                        let s = vertices.len() as u32;
                        indices.push(Point3::new(s + 0, s + 2, s + 1));
                        indices.push(Point3::new(s + 0, s + 3, s + 2));
                        indices.push(Point3::new(s + 4, s + 5, s + 6));
                        indices.push(Point3::new(s + 4, s + 6, s + 7));
                        indices.push(Point3::new(s + 7, s + 6, s + 2));
                        indices.push(Point3::new(s + 7, s + 2, s + 3));
                        indices.push(Point3::new(s + 4, s + 1, s + 5));
                        indices.push(Point3::new(s + 4, s + 0, s + 1));
                        indices.push(Point3::new(s + 6, s + 5, s + 1));
                        indices.push(Point3::new(s + 6, s + 1, s + 2));
                        indices.push(Point3::new(s + 7, s + 0, s + 4));
                        indices.push(Point3::new(s + 7, s + 3, s + 0));
                    }
                }
            }
        }

        (vertices, indices)
    }
}

impl Into<VoxelSet> for VoxelizedVolume {
    fn into(self) -> VoxelSet {
        let mut vset = VoxelSet::new();
        vset.origin = self.origin;
        vset.voxels
            .reserve((self.num_voxels_inside_surface + self.num_voxels_on_surface) as usize);
        vset.scale = self.scale;

        for i in 0..self.resolution.x {
            for j in 0..self.resolution.y {
                for k in 0..self.resolution.z {
                    let value = self.voxel(i, j, k);

                    if value == VoxelValue::PrimitiveInsideSurface {
                        let voxel = Voxel {
                            coords: Point3::new(i, j, k),
                            is_on_surface: false,
                        };
                        vset.voxels.push(voxel);
                    } else if value == VoxelValue::PrimitiveOnSurface {
                        let voxel = Voxel {
                            coords: Point3::new(i, j, k),
                            is_on_surface: true,
                        };
                        vset.voxels.push(voxel);
                    }
                }
            }
        }

        vset
    }
}

/*
fn traceRay(
    mesh: &RaycastMesh,
    start: Real,
    dir: &Vector3<Real>,
    inside_count: &mut u32,
    outside_count: &mut u32,
) {
    let out_t;
    let u;
    let v;
    let w;
    let face_sign;
    let face_index;
    let hit = raycast_mesh.raycast(start, dir, out_t, u, v, w, face_sign, face_index);

    if hit {
        if face_sign >= 0 {
            *inside_count += 1;
        } else {
            *outside_count += 1;
        }
    }
}


fn raycast_fill(volume: &Volume, raycast_mesh: &RaycastMesh) {
if !raycast_mesh {
    return;
}

let scale = volume.scale;
let bmin = volume.min_bb;

let i0 = volume.resolution[0];
let j0 = volume.resolution[1];
let k0 = volume.resolution[2];

for i in 0..i0 {
    for j in 0..j0 {
        for k in 0..k0 {
            let voxel = volume.get_voxel(i, j, k);

            if voxel != VoxelValue::PrimitiveOnSurface {
                let start = Vector3::new(
                    i as Real * scale + bmin[0],
                    j as Real * scale + bmin[1],
                    k as Real * scale + bmin[2],
                );

                let mut inside_count = 0;
                let mut outside_count = 0;

                let directions = [
                    Vector3::x(),
                    -Vector3::x(),
                    Vector3::y(),
                    -Vector3::y(),
                    Vector3::z(),
                    -Vector3::z(),
                ];

                for r in 0..6 {
                    traceRay(
                        raycast_mesh,
                        start,
                        &directions[r * 3],
                        &mut inside_count,
                        &mut outside_count,
                    );

                    // Early out if we hit the outside of the mesh
                    if outside_count != 0 {
                        break;
                    }

                    // Early out if we accumulated 3 inside hits
                    if inside_count >= 3 {
                        break;
                    }
                }

                if outside_count == 0 && inside_count >= 3 {
                    volume.set_voxel(i, j, k, VoxelValue::PrimitiveInsideSurface);
                } else {
                    volume.set_voxel(i, j, k, VoxelValue::PrimitiveOutsideSurface);
                }
            }
        }
    }
}
}
 */
