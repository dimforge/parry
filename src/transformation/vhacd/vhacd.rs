// Rust port, with modifications, of https://github.com/kmammou/v-hacd/blob/master/src/VHACD_Lib/src/VHACD.cpp
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

use crate::math::{Point, Real, Vector, DIM};
use crate::transformation::vhacd::VHACDParameters;
use crate::transformation::voxelization::{VoxelSet, VoxelizedVolume};
use std::sync::Arc;

#[cfg(feature = "dim2")]
type ConvexHull = Vec<Point<Real>>;
#[cfg(feature = "dim3")]
type ConvexHull = (Vec<Point<Real>>, Vec<[u32; 3]>);

#[derive(Copy, Clone, Debug)]
pub(crate) struct CutPlane {
    pub abc: Vector<Real>,
    pub d: Real,
    pub axis: u8,
    pub index: u32,
}

/// Approximate convex decomposition using the VHACD algorithm.
pub struct VHACD {
    // raycast_mesh: Option<RaycastMesh>,
    voxel_parts: Vec<VoxelSet>,
    volume_ch0: Real,
    max_concavity: Real,
}

impl VHACD {
    /// Decompose the given polyline (in 2D) or triangle mesh (in 3D).
    ///
    /// # Parameters
    /// * `params` - The parameters for the VHACD algorithm execution.
    /// * `points` - The vertex buffer of the polyline (in 2D) or triangle mesh (in 3D).
    /// * `indices` - The index buffer of the polyline (in 2D) or triangle mesh (in 3D).
    /// * `keep_voxel_to_primitives_map` - If set to `true` then a map between the voxels
    ///   computed during the decomposition, and the primitives (triangle or segment) they
    ///   intersect will be computed. This is required in order to compute the convex-hulls
    ///   using the original polyline/trimesh primitives afterwards (otherwise the convex
    ///   hulls resulting from the convex decomposition will use the voxels vertices).
    pub fn decompose(
        params: &VHACDParameters,
        points: &[Point<Real>],
        indices: &[[u32; DIM]],
        keep_voxel_to_primitives_map: bool,
    ) -> Self {
        // if params.project_hull_vertices || params.fill_mode == FillMode::RAYCAST_FILL {
        //     self.raycast_mesh =
        //         RaycastMesh::create_raycast_mesh(num_points, points, num_triangles, triangles);
        // }

        let voxelized = VoxelizedVolume::voxelize(
            points,
            indices,
            params.resolution,
            params.fill_mode,
            keep_voxel_to_primitives_map,
            // &self.raycast_mesh,
        );

        let mut result = Self::from_voxels(params, voxelized.into());

        let primitive_classes = Arc::new(result.classify_primitives(indices.len()));
        for part in &mut result.voxel_parts {
            part.primitive_classes = primitive_classes.clone();
        }

        result
    }

    /// Perform an approximate convex decomposition of a set of voxels.
    pub fn from_voxels(params: &VHACDParameters, voxels: VoxelSet) -> Self {
        let mut result = Self {
            // raycast_mesh: None,
            voxel_parts: Vec::new(),
            volume_ch0: 0.0,
            max_concavity: -Real::MAX,
        };

        result.do_compute_acd(params, voxels);
        result
    }

    /// The almost-convex voxelized parts computed by the VHACD algorithm.
    pub fn voxel_parts(&self) -> &[VoxelSet] {
        &self.voxel_parts
    }

    #[cfg(feature = "dim2")]
    fn compute_preferred_cutting_direction(eigenvalues: &Vector<Real>) -> (Vector<Real>, Real) {
        let vx = eigenvalues.y * eigenvalues.y;
        let vy = eigenvalues.x * eigenvalues.x;

        if vx < vy {
            let e = eigenvalues.y * eigenvalues.y;
            let dir = Vector::x();

            if e == 0.0 {
                (dir, 0.0)
            } else {
                (dir, 1.0 - vx / e)
            }
        } else {
            let e = eigenvalues.x * eigenvalues.x;
            let dir = Vector::y();

            if e == 0.0 {
                (dir, 0.0)
            } else {
                (dir, 1.0 - vy / e)
            }
        }
    }

    #[cfg(feature = "dim3")]
    fn compute_preferred_cutting_direction(eigenvalues: &Vector<Real>) -> (Vector<Real>, Real) {
        let vx = (eigenvalues.y - eigenvalues.z) * (eigenvalues.y - eigenvalues.z);
        let vy = (eigenvalues.x - eigenvalues.z) * (eigenvalues.x - eigenvalues.z);
        let vz = (eigenvalues.x - eigenvalues.y) * (eigenvalues.x - eigenvalues.y);

        if vx < vy && vx < vz {
            let e = eigenvalues.y * eigenvalues.y + eigenvalues.z * eigenvalues.z;
            let dir = Vector::x();

            if e == 0.0 {
                (dir, 0.0)
            } else {
                (dir, 1.0 - vx / e)
            }
        } else if vy < vx && vy < vz {
            let e = eigenvalues.x * eigenvalues.x + eigenvalues.z * eigenvalues.z;
            let dir = Vector::y();

            if e == 0.0 {
                (dir, 0.0)
            } else {
                (dir, 1.0 - vy / e)
            }
        } else {
            let e = eigenvalues.x * eigenvalues.x + eigenvalues.y * eigenvalues.y;
            let dir = Vector::z();

            if e == 0.0 {
                (dir, 0.0)
            } else {
                (dir, 1.0 - vz / e)
            }
        }
    }

    // TODO: this should be a method of VoxelSet.
    fn compute_axes_aligned_clipping_planes(
        vset: &VoxelSet,
        downsampling: u32,
        planes: &mut Vec<CutPlane>,
    ) {
        let min_v = vset.min_bb_voxels();
        let max_v = vset.max_bb_voxels();

        for dim in 0..DIM {
            let i0 = min_v[dim];
            let i1 = max_v[dim];

            for i in (i0..=i1).step_by(downsampling as usize) {
                let plane = CutPlane {
                    abc: Vector::ith(dim, 1.0),
                    axis: dim as u8,
                    d: -(vset.origin[dim] + (i as Real + 0.5) * vset.scale),
                    index: i,
                };

                planes.push(plane);
            }
        }
    }

    fn refine_axes_aligned_clipping_planes(
        vset: &VoxelSet,
        best_plane: &CutPlane,
        downsampling: u32,
        planes: &mut Vec<CutPlane>,
    ) {
        let min_v = vset.min_bb_voxels();
        let max_v = vset.max_bb_voxels();

        let best_id = best_plane.axis as usize;
        let i0 = min_v[best_id].max(best_plane.index.saturating_sub(downsampling));
        let i1 = max_v[best_id].min(best_plane.index + downsampling);

        for i in i0..=i1 {
            let plane = CutPlane {
                abc: Vector::ith(best_id, 1.0),
                axis: best_plane.axis,
                d: -(vset.origin[best_id] + (i as Real + 0.5) * vset.scale),
                index: i,
            };
            planes.push(plane);
        }
    }

    // Returns the best plane, and the min concavity.
    fn compute_best_clipping_plane(
        &self,
        input_voxels: &VoxelSet,
        input_voxels_ch: &ConvexHull,
        planes: &[CutPlane],
        preferred_cutting_direction: &Vector<Real>,
        w: Real,
        alpha: Real,
        beta: Real,
        convex_hull_downsampling: u32,
        params: &VHACDParameters,
    ) -> (CutPlane, Real) {
        let mut best_plane = planes[0];
        let mut min_concavity = Real::MAX;
        let mut i_best = -1;
        let mut min_total = Real::MAX;

        let mut left_ch;
        let mut right_ch;
        let mut left_ch_pts = Vec::new();
        let mut right_ch_pts = Vec::new();
        let mut left_voxels = VoxelSet::new();
        let mut right_voxels = VoxelSet::new();
        let mut on_surface_voxels = VoxelSet::new();

        input_voxels.select_on_surface(&mut on_surface_voxels);

        for (x, plane) in planes.iter().enumerate() {
            // Compute convex hulls.
            if params.convex_hull_approximation {
                right_ch_pts.clear();
                left_ch_pts.clear();

                on_surface_voxels.intersect(
                    plane,
                    &mut right_ch_pts,
                    &mut left_ch_pts,
                    convex_hull_downsampling * 32,
                );

                clip_mesh(
                    #[cfg(feature = "dim2")]
                    &input_voxels_ch,
                    #[cfg(feature = "dim3")]
                    &input_voxels_ch.0,
                    plane,
                    &mut right_ch_pts,
                    &mut left_ch_pts,
                );
                right_ch = convex_hull(&right_ch_pts);
                left_ch = convex_hull(&left_ch_pts);
            } else {
                on_surface_voxels.clip(plane, &mut right_voxels, &mut left_voxels);
                right_ch = right_voxels.compute_convex_hull(convex_hull_downsampling);
                left_ch = left_voxels.compute_convex_hull(convex_hull_downsampling);
            }

            let volume_left_ch = compute_volume(&left_ch);
            let volume_right_ch = compute_volume(&right_ch);

            // compute clipped volumes
            let (volume_left, volume_right) = input_voxels.compute_clipped_volumes(plane);
            let concavity_left = compute_concavity(volume_left, volume_left_ch, self.volume_ch0);
            let concavity_right = compute_concavity(volume_right, volume_right_ch, self.volume_ch0);
            let concavity = concavity_left + concavity_right;

            // compute cost
            let balance = alpha * (volume_left - volume_right).abs() / self.volume_ch0;
            let d = w * plane.abc.dot(preferred_cutting_direction);
            let symmetry = beta * d;
            let total = concavity + balance + symmetry;

            if total < min_total || (total == min_total && (x as i32) < i_best) {
                min_concavity = concavity;
                best_plane = *plane;
                min_total = total;
                i_best = x as i32;
            }
        }

        (best_plane, min_concavity)
    }

    fn process_primitive_set(
        &mut self,
        params: &VHACDParameters,
        first_iteration: bool,
        parts: &mut Vec<VoxelSet>,
        temp: &mut Vec<VoxelSet>,
        mut voxels: VoxelSet,
    ) {
        let volume = voxels.compute_volume(); // Compute the volume for this primitive set
        voxels.compute_bb(); // Compute the bounding box for this primitive set.
        let voxels_convex_hull = voxels.compute_convex_hull(params.convex_hull_downsampling); // Generate the convex hull for this primitive set.

        // Compute the volume of the convex hull
        let volume_ch = compute_volume(&voxels_convex_hull);

        // If this is the first iteration, store the volume of the base
        if first_iteration {
            self.volume_ch0 = volume_ch;
        }

        // Compute the concavity of this volume
        let concavity = compute_concavity(volume, volume_ch, self.volume_ch0);

        // Compute the volume error.
        if concavity > params.concavity {
            let eigenvalues = voxels.compute_principal_axes();
            let (preferred_cutting_direction, w) =
                Self::compute_preferred_cutting_direction(&eigenvalues);

            let mut planes = Vec::new();
            Self::compute_axes_aligned_clipping_planes(
                &voxels,
                params.plane_downsampling,
                &mut planes,
            );

            let (mut best_plane, mut min_concavity) = self.compute_best_clipping_plane(
                &voxels,
                &voxels_convex_hull,
                &planes,
                &preferred_cutting_direction,
                w,
                concavity * params.alpha,
                concavity * params.beta,
                params.convex_hull_downsampling,
                params,
            );

            if params.plane_downsampling > 1 || params.convex_hull_downsampling > 1 {
                let mut planes_ref = Vec::new();

                Self::refine_axes_aligned_clipping_planes(
                    &voxels,
                    &best_plane,
                    params.plane_downsampling,
                    &mut planes_ref,
                );

                let best = self.compute_best_clipping_plane(
                    &voxels,
                    &voxels_convex_hull,
                    &planes_ref,
                    &preferred_cutting_direction,
                    w,
                    concavity * params.alpha,
                    concavity * params.beta,
                    1, // convex_hull_downsampling = 1
                    params,
                );

                best_plane = best.0;
                min_concavity = best.1;
            }

            if min_concavity > self.max_concavity {
                self.max_concavity = min_concavity;
            }

            let mut best_left = VoxelSet::new();
            let mut best_right = VoxelSet::new();

            voxels.clip(&best_plane, &mut best_right, &mut best_left);

            temp.push(best_left);
            temp.push(best_right);
        } else {
            parts.push(voxels);
        }
    }

    fn do_compute_acd(&mut self, params: &VHACDParameters, mut voxels: VoxelSet) {
        let intersections = voxels.intersections.clone();
        let mut input_parts = Vec::new();
        let mut parts = Vec::new();
        let mut temp = Vec::new();
        input_parts.push(std::mem::replace(&mut voxels, VoxelSet::new()));

        let mut first_iteration = true;
        self.volume_ch0 = 1.0;

        // Compute the decomposition depth based on the number of convex hulls being requested.
        let mut hull_count = 2;
        let mut depth = 1;

        while params.max_convex_hulls > hull_count {
            depth += 1;
            hull_count *= 2;
        }

        // We must always increment the decomposition depth one higher than the maximum number of hulls requested.
        // The reason for this is as follows.
        // Say, for example, the user requests 32 convex hulls exactly.  This would be a decomposition depth of 5.
        // However, when we do that, we do *not* necessarily get 32 hulls as a result.  This is because, during
        // the recursive descent of the binary tree, one or more of the leaf nodes may have no concavity and
        // will not be split.  So, in this way, even with a decomposition depth of 5, you can produce fewer than
        // 32 hulls.  So, in this case, we would set the decomposition depth to 6 (producing up to as high as 64 convex
        // hulls). Then, the merge step which combines over-described hulls down to the user requested amount, we will end
        // up getting exactly 32 convex hulls as a result. We could just allow the artist to directly control the
        // decomposition depth directly, but this would be a bit too complex and the preference is simply to let them
        // specify how many hulls they want and derive the solution from that.
        depth += 1;

        for _ in 0..depth {
            if input_parts.is_empty() {
                break;
            }

            for input_part in input_parts.drain(..) {
                self.process_primitive_set(
                    params,
                    first_iteration,
                    &mut parts,
                    &mut temp,
                    input_part,
                );
                first_iteration = false;
            }

            std::mem::swap(&mut input_parts, &mut temp);
            // Note that temp is already clear because our previous for
            // loop used `drain`. However we call `clear` here explicitly
            // to make sure it still works if we remove the `drain` in the
            // future.
            temp.clear();
        }

        parts.append(&mut input_parts);
        self.voxel_parts = parts;

        for part in &mut self.voxel_parts {
            part.intersections = intersections.clone();
        }
    }

    // Returns a vector such that `result[i]` gives the index of the the voxelized convex part that
    // intersects it.
    //
    // If multiple convex parts intersect the same primitive, then `result[i]` is set to `u32::MAX`.
    // This is used to avoid some useless triangle/segment cutting when computing the exact convex hull
    // of a voxelized convex part.
    fn classify_primitives(&self, num_primitives: usize) -> Vec<u32> {
        if num_primitives == 0 {
            return Vec::new();
        }

        const NO_CLASS: u32 = u32::MAX - 1;
        const MULTICLASS: u32 = u32::MAX;

        let mut primitive_classes = Vec::new();
        primitive_classes.resize(num_primitives, NO_CLASS);

        for (ipart, part) in self.voxel_parts.iter().enumerate() {
            for voxel in &part.voxels {
                let range = voxel.intersections_range.0..voxel.intersections_range.1;
                for inter in &part.intersections[range] {
                    let class = &mut primitive_classes[*inter as usize];
                    if *class == NO_CLASS {
                        *class = ipart as u32;
                    } else if *class != ipart as u32 {
                        *class = MULTICLASS;
                    }
                }
            }
        }

        primitive_classes
    }

    /// Compute the intersections between the voxelized convex part of this decomposition,
    /// and all the primitives from the original decomposed polyline/trimesh,
    ///
    /// This will panic if `keep_voxel_to_primitives_map` was set to `false` when initializing
    /// `self`.
    pub fn compute_primitive_intersections(
        &self,
        points: &[Point<Real>],
        indices: &[[u32; DIM]],
    ) -> Vec<Vec<Point<Real>>> {
        self.voxel_parts
            .iter()
            .map(|part| part.compute_primitive_intersections(points, indices))
            .collect()
    }

    /// Compute the convex-hulls of the parts computed by this approximate convex-decomposition,
    /// taking into account the primitives from the original polyline/trimesh being decomposed.
    ///
    /// This will panic if `keep_voxel_to_primitives_map` was set to `false` when initializing
    /// `self`.
    #[cfg(feature = "dim2")]
    pub fn compute_exact_convex_hulls(
        &self,
        points: &[Point<Real>],
        indices: &[[u32; DIM]],
    ) -> Vec<Vec<Point<Real>>> {
        self.voxel_parts
            .iter()
            .map(|part| part.compute_exact_convex_hull(points, indices))
            .collect()
    }

    /// Compute the convex-hulls of the parts computed by this approximate convex-decomposition,
    /// taking into account the primitives from the original polyline/trimesh being decomposed.
    ///
    /// This will panic if `keep_voxel_to_primitives_map` was set to `false` when initializing
    /// `self`.
    #[cfg(feature = "dim3")]
    pub fn compute_exact_convex_hulls(
        &self,
        points: &[Point<Real>],
        indices: &[[u32; DIM]],
    ) -> Vec<(Vec<Point<Real>>, Vec<[u32; DIM]>)> {
        self.voxel_parts
            .iter()
            .map(|part| part.compute_exact_convex_hull(points, indices))
            .collect()
    }

    /// Compute the convex hulls of the voxelized approximately-convex parts
    /// computed by `self` on the voxelized model.
    ///
    /// Use `compute_exact_convex_hulls` instead if the original polyline/trimesh geometry
    /// needs to be taken into account.
    #[cfg(feature = "dim2")]
    pub fn compute_convex_hulls(&self, downsampling: u32) -> Vec<Vec<Point<Real>>> {
        let downsampling = downsampling.max(1);
        self.voxel_parts
            .iter()
            .map(|part| part.compute_convex_hull(downsampling))
            .collect()
    }

    /// Compute the convex hulls of the voxelized approximately-convex parts
    /// computed by `self` on the voxelized model.
    ///
    /// Use `compute_exact_convex_hulls` instead if the original polyline/trimesh geometry
    /// needs to be taken into account.
    #[cfg(feature = "dim3")]
    pub fn compute_convex_hulls(
        &self,
        downsampling: u32,
    ) -> Vec<(Vec<Point<Real>>, Vec<[u32; DIM]>)> {
        let downsampling = downsampling.max(1);
        self.voxel_parts
            .iter()
            .map(|part| part.compute_convex_hull(downsampling))
            .collect()
    }
}

fn compute_concavity(volume: Real, volume_ch: Real, volume0: Real) -> Real {
    (volume_ch - volume).abs() / volume0
}

fn clip_mesh(
    points: &[Point<Real>],
    plane: &CutPlane,
    positive_part: &mut Vec<Point<Real>>,
    negative_part: &mut Vec<Point<Real>>,
) {
    for pt in points {
        let d = plane.abc.dot(&pt.coords) + plane.d;

        if d > 0.0 {
            positive_part.push(*pt);
        } else if d < 0.0 {
            negative_part.push(*pt);
        } else {
            positive_part.push(*pt);
            negative_part.push(*pt);
        }
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

#[cfg(feature = "dim2")]
fn compute_volume(polygon: &[Point<Real>]) -> Real {
    if !polygon.is_empty() {
        crate::mass_properties::details::convex_polygon_area_and_center_of_mass(&polygon).0
    } else {
        0.0
    }
}

#[cfg(feature = "dim3")]
fn compute_volume(mesh: &(Vec<Point<Real>>, Vec<[u32; DIM]>)) -> Real {
    if !mesh.0.is_empty() {
        crate::mass_properties::details::convex_mesh_volume_and_center_of_mass_unchecked(
            &mesh.0, &mesh.1,
        )
        .0
    } else {
        0.0
    }
}
