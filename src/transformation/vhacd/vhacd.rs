use super::{Axis, Mesh, Plane, Volume, VoxelSet};
use crate::na::Isometry;
use crate::transformation::vhacd::Parameters;
use crate::utils;
use na::{Point3, Vector3};

fn find_minimum_element(d: &[Real], m: &mut Real, begin: i32, end: i32) -> i32 {
    let mut idx = -1;
    let mut min = Real::MAX;

    for i in begin..end {
        if d[i as usize] < min {
            idx = i;
            min = d[i as usize];
        }
    }

    *m = min;
    return idx;
}

pub struct VHACD {
    // raycast_mesh: Option<RaycastMesh>,
    convex_hulls: Vec<Mesh>,
    volume_ch0: Real,
    dim: u32,
    voxel_resolution: Real, // 1
    volume: Volume,
    pset: VoxelSet,
    max_concavity: Real,
}

impl VHACD {
    pub fn compute(
        &mut self,
        points: &[Point3<Real>],
        triangles: &[Point3<u32>],
        params: &Parameters,
    ) {
        self.compute_acd(points, triangles, params)
    }

    fn init(&mut self) {
        // self.raycast_mesh = None;
        self.dim = 64;
        self.volume = Volume::new();
        self.volume_ch0 = 0.0;
        self.pset = VoxelSet::new();
        self.max_concavity = -Real::MAX;
    }

    fn compute_acd(
        &mut self,
        points: &[Point3<Real>],
        triangles: &[Point3<u32>],
        params: &Parameters,
    ) {
        self.init();

        // if params.project_hull_vertices || params.fill_mode == FillMode::RAYCAST_FILL {
        //     self.raycast_mesh =
        //         RaycastMesh::create_raycast_mesh(num_points, points, num_triangles, triangles);
        // }

        self.voxelize_mesh(points, triangles, params);
        self.compute_primitive_set(params);
        self.do_compute_acd(params);
        // self.merge_convex_hulls(params);
        // self.simplify_convex_hulls(params);
    }

    fn voxelize_mesh(
        &mut self,
        points: &[Point3<Real>],
        triangles: &[Point3<u32>],
        params: &Parameters,
    ) {
        // Default dimensions is the cube root of the resolution provided times the
        // default voxel dimension of 64
        let a = (params.resolution as Real).powf(0.33);
        self.dim = (a * 1.5) as u32;

        // Minimum voxel resolution is 32x32x32
        if self.dim < 32 {
            self.dim = 32;
        }

        self.volume = Volume::new();
        self.volume.voxelize(
            &Isometry::identity(),
            points,
            triangles,
            self.dim,
            params.fill_mode,
            // &self.raycast_mesh,
        );

        self.voxel_resolution = self.volume.scale();
    }

    fn compute_primitive_set(&mut self, params: &Parameters) {
        self.pset = VoxelSet::new();
        self.volume.convert_to_voxel_set(&mut self.pset);
        self.volume = Volume::new(); // Free the volume memory.
    }

    // TODO: this should just be a method of VoxelSet.
    fn compute_preferred_cutting_direction(tset: &VoxelSet) -> (Vector3<Real>, Real) {
        let eigv = tset.eigenvalues();

        let vx = (eigv.y - eigv.z) * (eigv.y - eigv.z);
        let vy = (eigv.x - eigv.z) * (eigv.x - eigv.z);
        let vz = (eigv.x - eigv.y) * (eigv.x - eigv.y);

        if vx < vy && vx < vz {
            let e = eigv.y * eigv.y + eigv.z * eigv.z;
            let dir = Vector3::x();

            if e == 0.0 {
                (dir, 0.0)
            } else {
                (dir, 1.0 - vx / e)
            }
        } else if vy < vx && vy < vz {
            let e = eigv.x * eigv.x + eigv.z * eigv.z;
            let dir = Vector3::y();

            if e == 0.0 {
                (dir, 0.0)
            } else {
                (dir, 1.0 - vy / e)
            }
        } else {
            let e = eigv.x * eigv.x + eigv.y * eigv.y;
            let dir = Vector3::z();

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
        planes: &mut Vec<Plane>,
    ) {
        let min_v = vset.min_bb_voxels();
        let max_v = vset.max_bb_voxels();

        for dim in 0..3 {
            let i0 = min_v[dim];
            let i1 = max_v[dim];

            for i in (i0..=i1).step_by(downsampling as usize) {
                let plane = Plane {
                    abc: Vector3::ith(dim, 1.0),
                    axis: Axis::try_from(dim).unwrap(),
                    d: -(i as Real + 0.5),
                    index: i,
                };

                planes.push(plane);
            }
        }
    }

    fn refine_axes_aligned_clipping_planes(
        vset: &VoxelSet,
        best_plane: &Plane,
        downsampling: u32,
        planes: &mut Vec<Plane>,
    ) {
        let min_v = vset.min_bb_voxels();
        let max_v = vset.max_bb_voxels();

        let best_id = best_plane.axis as usize;
        let i0 = min_v[best_id].max(best_plane.index - downsampling);
        let i1 = max_v[best_id].min(best_plane.index + downsampling);

        for i in i0..=i1 {
            let plane = Plane {
                abc: Vector3::ith(best_id, 1.0),
                axis: best_plane.axis,
                d: -(i as Real + 0.5),
                index: i,
            };
            planes.push(plane);
        }
    }

    // Returns the best plane, and the min concavity.
    fn compute_best_clipping_plane(
        &self,
        input_pset: &VoxelSet,
        input_pset_ch: Option<&Mesh>,
        planes: &[Plane],
        preferred_cutting_direction: &Vector3<Real>,
        w: Real,
        alpha: Real,
        beta: Real,
        convex_hull_downsampling: u32,
        params: &Parameters,
    ) -> (Plane, Real) {
        let mut best_plane = planes[0];
        let mut min_concavity = Real::MAX;
        let mut i_best = -1;
        let mut done = 0;
        let mut min_total = Real::MAX;
        let mut min_balance = Real::MAX;
        let mut min_symmetry = Real::MAX;

        let mut left_ch = Mesh::new();
        let mut right_ch = Mesh::new();
        let mut left_ch_pts = Vec::new();
        let mut right_ch_pts = Vec::new();
        let mut left_pset = VoxelSet::new();
        let mut right_pset = VoxelSet::new();
        let mut on_surface_pset = VoxelSet::new();

        input_pset.select_on_surface(&mut on_surface_pset);

        for (x, plane) in planes.iter().enumerate() {
            right_ch.points.clear();
            left_ch.points.clear();
            right_ch.triangles.clear();
            left_ch.triangles.clear();

            // Compute convex hulls.
            if params.convex_hull_approximation {
                right_ch_pts.clear();
                left_ch_pts.clear();

                on_surface_pset.intersect(
                    plane,
                    &mut right_ch_pts,
                    &mut left_ch_pts,
                    convex_hull_downsampling * 32,
                );

                input_pset_ch
                    .unwrap()
                    .clip(plane, &mut right_ch_pts, &mut left_ch_pts);
                right_ch.compute_convex_hull(&right_ch_pts);
                left_ch.compute_convex_hull(&left_ch_pts);
            } else {
                on_surface_pset.clip(plane, &mut right_pset, &mut left_pset);
                right_pset.compute_convex_hull(&mut right_ch, convex_hull_downsampling);
                left_pset.compute_convex_hull(&mut left_ch, convex_hull_downsampling);
            }

            let volume_left_ch = left_ch.compute_volume();
            let volume_right_ch = right_ch.compute_volume();

            // compute clipped volumes
            let (volume_left, volume_right) = input_pset.compute_clipped_volumes(plane);
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
                min_balance = balance;
                min_symmetry = symmetry;
                best_plane = *plane;
                min_total = total;
                i_best = x as i32;
            }

            done += 1;
        }

        (best_plane, min_concavity)
    }

    fn process_primitive_set(
        &mut self,
        params: &Parameters,
        first_iteration: bool,
        parts: &mut Vec<VoxelSet>,
        temp: &mut Vec<VoxelSet>,
        max_concavity: Real,
        mut pset: VoxelSet,
    ) {
        let volume = pset.compute_volume(); // Compute the volume for this primitive set
        pset.compute_bb(); // Compute the bounding box for this primitive set.
        pset.compute_principal_axes(); // Compute the principle axes.
        let mut pset_convex_hull = Mesh::new();
        pset.compute_convex_hull(&mut pset_convex_hull, params.convex_hull_downsampling); // Generate the convex hull for this primitive set.

        // Compute the volume of the convex hull
        let volume_ch = pset_convex_hull.compute_volume();

        // If this is the first iteration, store the volume of the base
        if first_iteration {
            self.volume_ch0 = volume_ch;
        }

        // Compute the concavity of this volume
        let concavity = compute_concavity(volume, volume_ch, self.volume_ch0);

        // Compute the volume error
        if concavity > params.concavity {
            let (preferred_cutting_direction, w) = Self::compute_preferred_cutting_direction(&pset);

            let mut planes = Vec::new();
            Self::compute_axes_aligned_clipping_planes(
                &pset,
                params.plane_downsampling,
                &mut planes,
            );

            let (mut best_plane, mut min_concavity) = self.compute_best_clipping_plane(
                &pset,
                Some(&pset_convex_hull),
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
                    &pset,
                    &best_plane,
                    params.plane_downsampling,
                    &mut planes_ref,
                );

                let best = self.compute_best_clipping_plane(
                    &pset,
                    Some(&pset_convex_hull),
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

            pset.clip(&best_plane, &mut best_right, &mut best_left);
            temp.push(best_left);
            temp.push(best_right);
        } else {
            parts.push(pset);
        }
    }

    fn do_compute_acd(&mut self, params: &Parameters) {
        let mut input_parts = Vec::new();
        let mut parts = Vec::new();
        let mut temp = Vec::new();
        input_parts.push(std::mem::replace(&mut self.pset, VoxelSet::new()));

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

            let max_concavity = 0.0;

            for input_part in input_parts.drain(..) {
                self.process_primitive_set(
                    params,
                    first_iteration,
                    &mut parts,
                    &mut temp,
                    max_concavity,
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

        self.convex_hulls.clear();

        for part in &parts {
            let mut mesh = Mesh::new();
            part.compute_convex_hull(&mut mesh, params.convex_hull_downsampling);
            self.convex_hulls.push(mesh);
        }

        parts.clear();
    }

    /*
    fn simplify_convex_hull(&self, ch: &mut mesh, nvertices: usize, min_volume: Real) {
        if nvertices <= 4 {
            return;
        }

        let ic_hull = ICHull::new();

        if self.raycast_mesh {
            // We project these points onto the original source mesh to increase precision
            // The voxelization process drops floating point precision so returned data points are not exactly lying on the
            // surface of the original source mesh.
            // The first step is we need to compute the bounding box of the mesh we are trying to build a convex hull for.
            // From this bounding box, we compute the length of the diagonal to get a relative size and center for point
            // projection
            let nPoints = ch.get_num_points();
            let input_points = ch.get_points_buffer();
            let mut bmin = input_points[0];
            let mut bmax = input_points[1];

            for i in 1..npoints {
                let p = &input_points[i];
                p.update_min_max(bmin, bmax);
            }

            let center; // Vector3
            let diagonal_length = center.get_center(bmin, bmax); // Get the center of the bounding box

            // This is the error threshold for determining if we should use the raycast result data point vs. the voxelized
            // result.
            let point_distance_threshold = self.voxel_resolution * 4; // can only snap if the point is within 2 times the voxel
            // grid resolution

            // If a new point is within 1/100th the diagonal length of the bounding volume we do not add it.  To do so would
            // create a thin sliver in the resulting convex hull
            let snap_distance_threshold = diagonal_length * 0.01;
            let snap_distance_threshold_squared = snap_distance_threshold * snap_distance_threshold;

            // Allocate buffer for projected vertices
            let mut output_points = vec![Vector3::zeros(), npoints];
            let outCount = 0;

            for i in 0..npoints {
                let input_point = &mut input_points[i];
                let output_point = &mut output_points[outCount];
                // Compute the direction vector from the center of this mesh to the vertex
                let mut dir = input_point - center;
                // Normalize the direction vector.
                dir.normalize_mut();
                // Multiply times the diagonal length of the mesh
                dir *= diagonal_length;
                // Add the center back in again to get the destination point
                dir += center;
                // By default the output point is equal to the input point
                output_point = input_point;
                let closest_point = [0.0; 3];

                // Find the closest point on the source mesh which is within the voxel resolution. If found
                // we 'snap' the convex hull input point to this surface mesh point.
                if self.raycast_mesh.get_closest_point_within_distance(&input_point[0], point_distance_threshold, closest_point) {
                    output_point[0] = closest_point[0];
                    output_point[1] = closest_point[1];
                    output_point[2] = closest_point[2];
                }

                // Ok, before we add this point, we do not want to create points which are extremely close to each other.
                // This will result in tiny sliver triangles which are really bad for collision detection.
                let mut found_nearby_point = false;

                for j in 0..out_count {
                    // If this new point is extremely close to an existing point, we do not add it!
                    let squared_distance = output_points[j].get_distance_squared(output_point);

                    if squared_distance < snap_distance_threshold_squared {
                        found_nearby_point = true;
                        break;
                    }
                }
                if !found_nearby_point {
                    outCount += 1;
                }
            }

            ic_hull.add_points(output_points, outCount);
        } else {
            ic_hull.add_points(ch.get_points_buffer(), ch.get_num_points());
        }

        ic_hull.process(nvertices, min_volume);

        let mesh = &mut ic_hull.GetMesh();
        let num_triangles = mesh.get_num_triangles();
        let num_vertices = mesh.get_num_vertices();
        ch.resize_points(num_vertices);
        ch.resize_triangles(num_triangles);
        mesh.get_ifs(ch.get_points_buffer(), ch.get_triangles_buffer());
    }


    fn merge_convex_hulls(&mut self, params: &Parameters) {
        // Get the current number of convex hulls
        let num_convex_hulls = self.convex_hulls.len();
        // Iteration counter
        let iteration = 0;

        // While we have more than at least one convex hull
        if num_convex_hulls > 1 {
            // Get the gamma error threshold for when to exit
            let mut pts = Vec::new();
            let mut combined_ch = Mesh::new();

            // Populate the cost matrix
            let mut idx = 0;
            let cost_matrix = Vec::new();
            cost_matrix.resize(((num_convex_hulls * num_convex_hulls) - num_convex_hulls) >> 1);
            for p1 in 1..num_convex_hulls {
                let volume1 = self.convex_hulls[p1].compute_volume();

                for p2 in 0..p1 {
                    compute_convex_hull(self.convex_hulls[p1], self.convex_hulls[p2], pts, &combined_ch);
                    cost_matrix[idx++] = compute_concavity(
                    volume1 + self.convex_hulls[p2].compute_volume(), combined_ch.compute_volume(), self.volume_ch0);
                }
            }

            // Until we cant merge below the maximum cost
            let cost_size = self.convex_hulls.len();

            loop {
                // Search for lowest cost
                let bestCost = Real::MAX;

                let addr = find_minimum_element(cost_matrix.Data(), &bestCost, 0, cost_matrix.len());
                if (cost_size - 1) < params.max_convex_hulls {
                    break;
                }

                let addrI = sqrt(1 + (8 * addr))) as i32 - 1 >> 1;
                let p1 = addrI + 1;
                let p2 = addr - ((addrI * (addrI + 1)) >> 1);

                assert!(p1 >= 0);
                assert!(p2 >= 0);
                assert!(p1 < cost_size);
                assert!(p2 < cost_size);

                // Make the lowest cost row and column into a new hull
                let mut cch = Mesh::new();

                compute_convex_hull(self.convex_hulls[p1], self.convex_hulls[p2], pts, cch);
                self.convex_hulls[p2] = cch;

                std::swap(self.convex_hulls[p1], self.convex_hulls[self.convex_hulls.len() - 1]);

                self.convex_hulls.pop(); // pop_back

                cost_size = cost_size - 1;

                // Calculate costs versus the new hull
                let row_idx = ((p2 - 1) * p2) >> 1;
                let volume1 = self.convex_hulls[p2].compute_volume();
                let mut i = 0;

                while i < p2 {
                    compute_convex_hull(self.convex_hulls[p2], self.convex_hulls[i], pts, &combined_ch);
                    cost_matrix[row_idx++] = compute_concavity(
                    volume1 + self.convex_hulls[i].compute_volume(), combined_ch.compute_volume(), self.volume_ch0);
                    i += 1;
                }

                row_idx += p2;
                i = p2 + 1;

                while i < cost_size {
                    compute_convex_hull(self.convex_hulls[p2], self.convex_hulls[i], pts, &combined_ch);
                    cost_matrix[row_idx] = compute_concavity(
                        volume1 + self.convex_hulls[i].compute_volume(),
                        combined_ch.compute_volume(),
                        self.volume_ch0
                    );
                    row_idx += i;
                    assert!(row_idx >= 0);
                    i += 1;
                }

                // Move the top column in to replace its space
                let erase_idx = ((cost_size - 1) * cost_size) >> 1;

                if p1 < cost_size {
                    row_idx = (addrI * p1) >> 1;
                    let top_row = erase_idx;

                    for i in 0..p1 {
                        if i != p2 {
                            cost_matrix[row_idx] = cost_matrix[top_row];
                        }

                        row_idx += 1;
                        top_row += 1;
                    }

                    top_row += 1;
                    row_idx += p1;

                    for i in p1 + 1..cost_size + 1 {
                        cost_matrix[row_idx] = cost_matrix[top_row];
                        top_row += 1;
                        row_idx += i;
                        assert!(row_idx >= 0);
                    }
                }

                cost_matrix.resize(erase_idx);
            }
        }
    }


    fn compute_center_of_mass(&self, center_of_mass: Point3<Real>) -> bool {
        let mut ret = false;

        center_of_mass[0] = 0;
        center_of_mass[1] = 0;
        center_of_mass[2] = 0;

        // Get number of convex hulls in the result
        let hull_count = get_num_convex_hulls();

        if hull_count {
            ret = true;
            let mut total_volume = 0;
            // Initialize the center of mass to zero
            *center_of_mass = Point3::origin();

            // Compute the total volume of all convex hulls
            for i in 0..hull_count {
                let mut ch = ConvexHull::new();
                get_convex_hull(i, &mut ch);
                total_volume += ch.volume;
            }

            // compute the reciprocal of the total volume
            let recip_volume = 1.0 / total_volume;

            // Add in the weighted by volume average of the center point of each convex hull
            for i in 0..hull_count {
                let mut ch = ConvexHull::new();
                get_convex_hull(i, &mut ch);
                let ratio = ch.volume * recip_volume;
                *center_of_mass += ch.center * ratio;
            }
        }

        ret
    }
     */
}

fn compute_local_concavity(volume: Real, volume_ch: Real) -> Real {
    compute_concavity(volume, volume_ch, volume_ch)
}

fn compute_concavity(volume: Real, volume_ch: Real, volume0: Real) -> Real {
    (volume_ch - volume).abs() / volume0
}

/*
struct PrimitiveSetBase<'a> {
    parameters: &'a Parameters,
    parts: Vec<&'a mut PrimitiveSet>,
    temp: Vec<&'a mut PrimitiveSet>,
    vhacd: &'a mut VHACD,
}

impl PrimitiveSetBase {
    fn get_concavity(&self) -> Real {
        self.parameters.concavity
    }

    fn get_alpha(&self) -> Real {
        self.parameters.alpha
    }

    fn get_beta(&self) -> Real {
        self.parameters.beta
    }

    fn get_convex_hull_downsampling(&self) -> u32 {
        self.parameters.convex_hull_downsampling
    }

    fn get_plane_downsampling(&self) -> u32 {
        self.parameters.plane_downsampling
    }

    fn refresh_concavity(&mut self, min_concavity: Real) {
        if *self.max_concavity < min_concavity {
            *self.max_concavity = min_concavity;
        }
    }

    fn push_temp(&mut self, pset: &mut PrimitiveSet) {
        self.temp.push(pset);
    }

    fn push_parts(&mut self, pset: &mut PrimitiveSet) {
        self.parts.push(pset);
    }
}


fn add_points(mesh: &Mesh, pts: &mut Vec<Vector3<Real>>) {
    const int32_t n = (int32_t)mesh.get_num_points();
    for i in 0..n {
        pts.push(mesh.get_point(i));
    }
}

fn compute_convex_hull(ch1: &Mesh, ch2: &Mesh, pts: &mut [Vector3<Real>], combined_ch: &mut Mesh) {
    pts.resize(0);

    add_points(ch1, pts);
    add_points(ch2, pts);

    btConvexHullComputer ch;
    ch.compute(pts, -1.0, -1.0);

    combined_ch.resize_points(0);
    combined_ch.resize_triangles(0);

    for v in 0..ch.vertices.len()  {
        combined_ch.add_point(Vector3::new(ch.vertices[v].x, ch.vertices[v].y, ch.vertices[v].z));
    }

    let nt = ch.faces.len();

    for t in 0..nt {
        let source_edge = &ch.edges[ch.faces[t]];
        let a = source_edge.get_source_vertex();
        let b = source_edge.get_target_vertex();

        let edge = source_edge.get_next_edge_of_face();

        let c = edge.get_target_vertex();

        while c != a {
            combined_ch.add_triangle(Vector3::new(a, b, c));
            edge = edge.get_next_edge_of_face();
            b = c;
            c = edge.get_target_vertex();
        }
    }
}
*/
