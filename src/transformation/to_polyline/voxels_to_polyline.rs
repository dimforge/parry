use crate::bounding_volume::Aabb;
use crate::math::{Point, Real, Vector};
use crate::shape::{VoxelType, Voxels};
use na::{self, Point2};

impl Voxels {
    /// Converts the outlines of this set of voxels as a polyline.
    pub fn to_polyline(&self) -> (Vec<Point2<Real>>, Vec<[u32; 2]>) {
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

    /// Outlines this voxels shape with segments.
    pub fn iter_outline(&self, mut f: impl FnMut(Point<Real>, Point<Real>)) {
        // TODO: move this as a new method: Voxels::to_outline?
        let radius = self.voxel_size() / 2.0;
        let aabb = Aabb::from_half_extents(Point::origin(), Vector::repeat(radius));
        let vtx = aabb.vertices();

        for (_, center, vox_data) in self.centers() {
            match vox_data.voxel_type() {
                VoxelType::Vertex => {
                    let mask = vox_data.feature_mask();

                    for edge in Aabb::FACES_VERTEX_IDS {
                        if mask & (1 << edge.0) != 0 || mask & (1 << edge.1) != 0 {
                            f(center + vtx[edge.0].coords, center + vtx[edge.1].coords);
                        }
                    }
                }
                VoxelType::Face => {
                    let vtx = aabb.vertices();
                    let mask = vox_data.feature_mask();

                    for (i, edge) in Aabb::FACES_VERTEX_IDS.iter().enumerate() {
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
