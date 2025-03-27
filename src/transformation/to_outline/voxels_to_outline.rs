use crate::bounding_volume::Aabb;
use crate::math::{Point, Real, Vector};
use crate::shape::{VoxelType, Voxels};
use crate::transformation::utils;
use na::{self, Point3, Vector3};

impl Voxels {
    pub fn to_outline(&self) -> (Vec<Point<Real>>, Vec<[u32; 2]>) {
        let mut points = vec![];
        self.iter_outline(|a, b| {
            points.push(a);
            points.push(b);
        });
        let indices = (0..points.len() as u32 / 2)
            .map(|i| [i * 2, i * 2 + 1])
            .collect();
        (points, indices)
    }

    /// Outlines this voxels shape using polylines.
    pub fn iter_outline(&self, mut f: impl FnMut(Point<Real>, Point<Real>)) {
        // TODO: move this as a new method: Voxels::to_outline?
        let radius = self.voxel_size() / 2.0;
        let aabb = Aabb::from_half_extents(Point::origin(), Vector::repeat(radius));
        let vtx = aabb.vertices();

        /*
        let aabb_oct = Aabb::from_half_extents(Point::origin(), Vector::repeat(radius / 2.0));
        let vtx_oct = aabb_oct.vertices();

        for (center, vox_data) in self.centers() {
            for k in 0..8 {
                let mask = (vox_data.octant_mask() >> (3 * k)) & 0b0111;
                match mask {
                    OctantPattern::EDGE_X | OctantPattern::EDGE_Y | OctantPattern::EDGE_Z => {
                        let i = mask - OctantPattern::EDGE_X;
                        let axis = Vector::ith(i as usize, 1.0);
                        f(
                            center + vtx_oct[k].coords + axis * radius / 2.0,
                            center + vtx_oct[k].coords - axis * radius / 2.0,
                        );
                    }
                    _ => {}
                }
            }

            for (i, edge) in Aabb::EDGES_VERTEX_IDS.iter().enumerate() {
                let mask0 = (vox_data.octant_mask() >> (3 * edge.0)) & 0b0111;
                let mask1 = (vox_data.octant_mask() >> (3 * edge.1)) & 0b0111;
                let mid = na::center(&vtx_oct[edge.0], &vtx_oct[edge.1]);

                if mask0 == OctantPattern::VERTEX {
                    f(center + vtx_oct[edge.0].coords, center + mid.coords);
                }

                if mask1 == OctantPattern::VERTEX {
                    f(center + mid.coords, center + vtx_oct[edge.1].coords);
                }
            }
        }
        */

        for (center, vox_data) in self.centers() {
            match vox_data.voxel_type() {
                VoxelType::Vertex => {
                    let mask = vox_data.feature_mask();

                    for edge in Aabb::EDGES_VERTEX_IDS {
                        if mask & (1 << edge.0) != 0 || mask & (1 << edge.1) != 0 {
                            f(center + vtx[edge.0].coords, center + vtx[edge.1].coords);
                        }
                    }
                }
                VoxelType::Edge => {
                    let vtx = aabb.vertices();
                    let mask = vox_data.feature_mask();

                    for (i, edge) in Aabb::EDGES_VERTEX_IDS.iter().enumerate() {
                        if mask & (1 << i) != 0 {
                            f(center + vtx[edge.0].coords, center + vtx[edge.1].coords);
                        }
                    }
                }
                _ => {}
            }
        }
    }
}
