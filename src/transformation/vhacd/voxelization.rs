use super::mesh::{Axis, Mesh, Plane};
use super::FillMode;
use super::VHACDParameters;
use crate::bounding_volume::AABB;
use crate::math::Real;
use crate::na::{Isometry3, SymmetricEigen};
use crate::query;
use crate::shape::Triangle;
use na::{Matrix3, Point3, Vector3};

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum VoxelValue {
    PrimitiveUndefined,
    PrimitiveOutsideSurfaceToWalk,
    PrimitiveOutsideSurface,
    PrimitiveInsideSurface,
    PrimitiveOnSurface,
}

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
    min_bb_voxels: Point3<u32>,
    max_bb_voxels: Point3<u32>,
    min_bb_pts: Point3<Real>,
    max_bb_pts: Point3<Real>,
    barycenter: Point3<u32>,
    barycenter_pca: Point3<Real>,
    pub scale: Real,
    unit_volume: Real,
    num_voxels_on_surface: u32,
    num_voxels_inside_surface: u32,
    voxels: Vec<Voxel>,
    eigenvalues: Vector3<Real>,
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

    pub fn compute_convex_hull(&self, mesh: &mut Mesh, sampling: u32) {
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
        mesh.compute_convex_hull(&points);
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

    fn compute_exterior_points(
        &self,
        plane: &Plane,
        mesh: &Mesh,
        exterior_pts: &mut Vec<Point3<Real>>,
    ) {
        let num_voxels = self.voxels.len();
        if num_voxels == 0 {
            return;
        }

        for v in 0..num_voxels {
            let voxel = self.voxels[v];
            let pt = self.get_voxel_point(&voxel);
            let d = plane.abc.dot(&pt.coords) + plane.d;

            if d >= 0.0 {
                if !mesh.is_inside(&pt) {
                    self.map_voxel_points(&voxel, |p| exterior_pts.push(p));
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
    fn convert(&self, mesh: &mut Mesh, value: VoxelValue) {
        for voxel in &self.voxels {
            if voxel.data == value {
                let s = mesh.points.len() as u32;
                self.map_voxel_points(voxel, |p| mesh.points.push(p));

                mesh.triangles.push(Point3::new(s + 0, s + 2, s + 1));
                mesh.triangles.push(Point3::new(s + 0, s + 3, s + 2));
                mesh.triangles.push(Point3::new(s + 4, s + 5, s + 6));
                mesh.triangles.push(Point3::new(s + 4, s + 6, s + 7));
                mesh.triangles.push(Point3::new(s + 7, s + 6, s + 2));
                mesh.triangles.push(Point3::new(s + 7, s + 2, s + 3));
                mesh.triangles.push(Point3::new(s + 4, s + 1, s + 5));
                mesh.triangles.push(Point3::new(s + 4, s + 0, s + 1));
                mesh.triangles.push(Point3::new(s + 6, s + 5, s + 1));
                mesh.triangles.push(Point3::new(s + 6, s + 1, s + 2));
                mesh.triangles.push(Point3::new(s + 7, s + 0, s + 4));
                mesh.triangles.push(Point3::new(s + 7, s + 3, s + 0));
            }
        }
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

pub struct Volume {
    dim: Point3<u32>,
    min_bb: Point3<Real>,
    max_bb: Point3<Real>,
    num_voxels_on_surface: u32,
    num_voxels_inside_surface: u32,
    num_voxels_outside_surface: u32,
    scale: Real,
    data: Vec<VoxelValue>,
}

impl Volume {
    pub fn new() -> Self {
        Volume {
            dim: Point3::origin(),
            min_bb: Point3::origin(),
            max_bb: Point3::new(1.0, 1.0, 1.0),
            num_voxels_on_surface: 0,
            num_voxels_inside_surface: 0,
            num_voxels_outside_surface: 0,
            scale: 1.0,
            data: Vec::new(),
        }
    }

    pub fn resolution(&self) -> Point3<u32> {
        self.dim
    }

    pub fn scale(&self) -> Real {
        self.scale
    }

    pub fn voxelize(
        &mut self,
        transform: &Isometry3<Real>,
        points: &[Point3<Real>],
        triangles: &[Point3<u32>],
        dim: u32,
        fill_mode: FillMode,
    ) {
        if points.is_empty() {
            return;
        }

        let aabb = crate::bounding_volume::point_cloud_aabb(transform, points);
        self.min_bb = aabb.mins;
        self.max_bb = aabb.maxs;

        let d = self.max_bb - self.min_bb;
        let r;

        if d[0] > d[1] && d[0] > d[2] {
            r = d[0];
            self.dim[0] = dim;
            self.dim[1] = 2 + (dim as Real * d[1] / d[0]) as u32;
            self.dim[2] = 2 + (dim as Real * d[2] / d[0]) as u32;
        } else if d[1] > d[0] && d[1] > d[2] {
            r = d[1];
            self.dim[1] = dim;
            self.dim[0] = 2 + (dim as Real * d[0] / d[1]) as u32;
            self.dim[2] = 2 + (dim as Real * d[2] / d[1]) as u32;
        } else {
            r = d[2];
            self.dim[2] = dim;
            self.dim[0] = 2 + (dim as Real * d[0] / d[2]) as u32;
            self.dim[1] = 2 + (dim as Real * d[1] / d[2]) as u32;
        }

        self.scale = r / (dim as Real - 1.0);
        let inv_scale = (dim as Real - 1.0) / r;
        self.allocate();

        self.num_voxels_on_surface = 0;
        self.num_voxels_on_surface = 0;
        self.num_voxels_inside_surface = 0;
        self.num_voxels_outside_surface = 0;

        let mut tri_pts = [Point3::origin(); 3];
        let box_half_size = Vector3::repeat(0.5);
        let mut i0 = 0;
        let mut i1 = 0;
        let mut j0 = 0;
        let mut j1 = 0;
        let mut k0 = 0;
        let mut k1 = 0;

        for tri in triangles {
            // Find the range of voxels potentially intersecting the triangle.
            for c in 0..3 {
                let pt = transform * points[tri[c] as usize];
                tri_pts[c] = (pt - self.min_bb.coords) * inv_scale;

                let i = (tri_pts[c].x + 0.5) as u32;
                let j = (tri_pts[c].y + 0.5) as u32;
                let k = (tri_pts[c].z + 0.5) as u32;

                assert!(i < self.dim[0] && j < self.dim[1] && k < self.dim[2]);

                if c == 0 {
                    i0 = i;
                    i1 = i;
                    j0 = j;
                    j1 = j;
                    k0 = k;
                    k1 = k;
                } else {
                    if i < i0 {
                        i0 = i;
                    }
                    if j < j0 {
                        j0 = j;
                    }
                    if k < k0 {
                        k0 = k;
                    }
                    if i > i1 {
                        i1 = i;
                    }
                    if j > j1 {
                        j1 = j;
                    }
                    if k > k1 {
                        k1 = k;
                    }
                }
            }

            if i0 > 0 {
                i0 -= 1;
            }

            if j0 > 0 {
                j0 -= 1;
            }

            if k0 > 0 {
                k0 -= 1;
            }

            if i1 < self.dim.x {
                i1 += 1;
            }

            if j1 < self.dim.y {
                j1 += 1;
            }

            if k1 < self.dim.z {
                k1 += 1;
            }

            // Determine exactly what voxel intersect the triangle.
            for i in i0..i1 {
                for j in j0..j1 {
                    for k in k0..k1 {
                        let value = self.voxel_mut(i, j, k);

                        if *value == VoxelValue::PrimitiveUndefined {
                            let triangle = Triangle::from(tri_pts);
                            let aabb = AABB::from_half_extents(
                                Point3::new(i as Real, j as Real, k as Real),
                                box_half_size,
                            );

                            if query::details::intersection_test_aabb_triangle(&aabb, &triangle) {
                                *value = VoxelValue::PrimitiveOnSurface;
                                self.num_voxels_on_surface += 1;
                            }
                        }
                    }
                }
            }
        }

        match fill_mode {
            FillMode::SurfaceOnly => {
                for value in &mut self.data {
                    if *value != VoxelValue::PrimitiveOnSurface {
                        *value = VoxelValue::PrimitiveOutsideSurface
                    }
                }
            }
            FillMode::FloodFill => {
                self.mark_outside_surface(0, 0, 0, self.dim[0], self.dim[1], 1);
                self.mark_outside_surface(
                    0,
                    0,
                    self.dim[2] - 1,
                    self.dim[0],
                    self.dim[1],
                    self.dim[2],
                );
                self.mark_outside_surface(0, 0, 0, self.dim[0], 1, self.dim[2]);
                self.mark_outside_surface(
                    0,
                    self.dim[1] - 1,
                    0,
                    self.dim[0],
                    self.dim[1],
                    self.dim[2],
                );
                self.mark_outside_surface(0, 0, 0, 1, self.dim[1], self.dim[2]);
                self.mark_outside_surface(
                    self.dim[0] - 1,
                    0,
                    0,
                    self.dim[0],
                    self.dim[1],
                    self.dim[2],
                );
                self.fill_outside_surface();
                self.fill_inside_surface();
            }
        }
    }

    fn allocate(&mut self) {
        let len = self.dim[0] * self.dim[1] * self.dim[2];
        self.data
            .resize(len as usize, VoxelValue::PrimitiveUndefined);
    }

    fn voxel_index(&self, i: u32, j: u32, k: u32) -> u32 {
        i + j * self.dim[0] + k * self.dim[0] * self.dim[1]
    }

    fn voxel_mut(&mut self, i: u32, j: u32, k: u32) -> &mut VoxelValue {
        let idx = self.voxel_index(i, j, k);
        &mut self.data[idx as usize]
    }

    pub fn voxel(&self, i: u32, j: u32, k: u32) -> VoxelValue {
        let idx = self.voxel_index(i, j, k);
        self.data[idx as usize]
    }

    fn free(&mut self) {
        self.data = Vec::new();
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
        let mut voxels_walked = 0;
        let i0 = self.dim[0];
        let j0 = self.dim[1];
        let k0 = self.dim[2];

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
        for i in 0..self.dim.x {
            for j in 0..self.dim.y {
                for k in 0..self.dim.z {
                    let v = self.voxel_mut(i, j, k);
                    if *v == VoxelValue::PrimitiveUndefined {
                        *v = VoxelValue::PrimitiveInsideSurface;
                        self.num_voxels_inside_surface += 1;
                    }
                }
            }
        }
    }

    pub fn convert_to_mesh(&self, mesh: &mut Mesh, value: VoxelValue) {
        let i0 = self.dim[0];
        let j0 = self.dim[1];
        let k0 = self.dim[2];

        for i in 0..i0 {
            for j in 0..j0 {
                for k in 0..k0 {
                    let voxel = self.voxel(i, j, k);

                    if voxel == value {
                        let i = i as Real;
                        let j = j as Real;
                        let k = k as Real;

                        let p0 = Vector3::new(
                            (i - 0.5) * self.scale,
                            (j - 0.5) * self.scale,
                            (k - 0.5) * self.scale,
                        );
                        let p1 = Vector3::new(
                            (i + 0.5) * self.scale,
                            (j - 0.5) * self.scale,
                            (k - 0.5) * self.scale,
                        );
                        let p2 = Vector3::new(
                            (i + 0.5) * self.scale,
                            (j + 0.5) * self.scale,
                            (k - 0.5) * self.scale,
                        );
                        let p3 = Vector3::new(
                            (i - 0.5) * self.scale,
                            (j + 0.5) * self.scale,
                            (k - 0.5) * self.scale,
                        );
                        let p4 = Vector3::new(
                            (i - 0.5) * self.scale,
                            (j - 0.5) * self.scale,
                            (k + 0.5) * self.scale,
                        );
                        let p5 = Vector3::new(
                            (i + 0.5) * self.scale,
                            (j - 0.5) * self.scale,
                            (k + 0.5) * self.scale,
                        );
                        let p6 = Vector3::new(
                            (i + 0.5) * self.scale,
                            (j + 0.5) * self.scale,
                            (k + 0.5) * self.scale,
                        );
                        let p7 = Vector3::new(
                            (i - 0.5) * self.scale,
                            (j + 0.5) * self.scale,
                            (k + 0.5) * self.scale,
                        );

                        let s = mesh.points.len() as u32;

                        mesh.points.push(self.min_bb + p0);
                        mesh.points.push(self.min_bb + p1);
                        mesh.points.push(self.min_bb + p2);
                        mesh.points.push(self.min_bb + p3);
                        mesh.points.push(self.min_bb + p4);
                        mesh.points.push(self.min_bb + p5);
                        mesh.points.push(self.min_bb + p6);
                        mesh.points.push(self.min_bb + p7);

                        mesh.triangles.push(Point3::new(s + 0, s + 2, s + 1));
                        mesh.triangles.push(Point3::new(s + 0, s + 3, s + 2));
                        mesh.triangles.push(Point3::new(s + 4, s + 5, s + 6));
                        mesh.triangles.push(Point3::new(s + 4, s + 6, s + 7));
                        mesh.triangles.push(Point3::new(s + 7, s + 6, s + 2));
                        mesh.triangles.push(Point3::new(s + 7, s + 2, s + 3));
                        mesh.triangles.push(Point3::new(s + 4, s + 1, s + 5));
                        mesh.triangles.push(Point3::new(s + 4, s + 0, s + 1));
                        mesh.triangles.push(Point3::new(s + 6, s + 5, s + 1));
                        mesh.triangles.push(Point3::new(s + 6, s + 1, s + 2));
                        mesh.triangles.push(Point3::new(s + 7, s + 0, s + 4));
                        mesh.triangles.push(Point3::new(s + 7, s + 3, s + 0));
                    }
                }
            }
        }
    }

    pub fn convert_to_voxel_set(&self, vset: &mut VoxelSet) {
        vset.min_bb = self.min_bb;
        vset.voxels
            .reserve((self.num_voxels_inside_surface + self.num_voxels_on_surface) as usize);
        vset.scale = self.scale;
        vset.unit_volume = self.scale * self.scale * self.scale;
        vset.num_voxels_on_surface = 0;
        vset.num_voxels_inside_surface = 0;

        for i in 0..self.dim.x {
            for j in 0..self.dim.y {
                for k in 0..self.dim.z {
                    let value = self.voxel(i, j, k);

                    if value == VoxelValue::PrimitiveInsideSurface {
                        let voxel = Voxel {
                            coords: Point3::new(i, j, k),
                            data: VoxelValue::PrimitiveInsideSurface,
                        };
                        vset.voxels.push(voxel);
                        vset.num_voxels_inside_surface += 1;
                    } else if value == VoxelValue::PrimitiveOnSurface {
                        let voxel = Voxel {
                            coords: Point3::new(i, j, k),
                            data: VoxelValue::PrimitiveOnSurface,
                        };
                        vset.voxels.push(voxel);
                        vset.num_voxels_on_surface += 1;
                    }
                }
            }
        }
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

let i0 = volume.dim[0];
let j0 = volume.dim[1];
let k0 = volume.dim[2];

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
