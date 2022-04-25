use crate::bounding_volume::AABB;
use crate::math::{Point, Real, Vector};
use crate::shape::RoundCuboid;
use crate::transformation::utils;

impl RoundCuboid {
    /// Discretize the boundary of this cuboid as a polyline.
    pub fn to_outline(&self, nsubdivs: u32) -> (Vec<Point<Real>>, Vec<[u32; 2]>) {
        let aabb = AABB::from_half_extents(Point::origin(), self.inner_shape.half_extents);
        let vtx = aabb.vertices();
        let fidx = AABB::FACES_VERTEX_IDS;
        let x = Vector::x() * self.border_radius;
        let y = Vector::y() * self.border_radius;
        let z = Vector::z() * self.border_radius;

        #[rustfmt::skip]
        let vtx = vec![
            vtx[fidx[0].0] + x, vtx[fidx[0].1] + x, vtx[fidx[0].2] + x, vtx[fidx[0].3] + x, vtx[fidx[0].4] + x,
            vtx[fidx[1].0] - x, vtx[fidx[1].1] - x, vtx[fidx[1].2] - x, vtx[fidx[1].3] - x, vtx[fidx[1].4] - x,
            vtx[fidx[2].0] + y, vtx[fidx[2].1] + y, vtx[fidx[2].2] + y, vtx[fidx[2].3] + y, vtx[fidx[2].4] + y,
            vtx[fidx[3].0] - y, vtx[fidx[3].1] - y, vtx[fidx[3].2] - y, vtx[fidx[3].3] - y, vtx[fidx[3].4] - y,
            vtx[fidx[4].0] + z, vtx[fidx[4].1] + z, vtx[fidx[4].2] + z, vtx[fidx[4].3] + z, vtx[fidx[4].4] + z,
            vtx[fidx[5].0] - z, vtx[fidx[5].1] - z, vtx[fidx[5].2] - z, vtx[fidx[5].3] - z, vtx[fidx[5].4] - z,
        ];

        let mut idx = vec![];
        for i in 0..6 {
            idx.push([i * 4 + 0, i * 4 + 1]);
            idx.push([i * 4 + 1, i * 4 + 2]);
            idx.push([i * 4 + 2, i * 4 + 3]);
            idx.push([i * 4 + 3, i * 4 + 0]);
        }

        let arcs = [
            (0, [4, 11, 20]),
            (0, [0, 12, 21]),
            (0, [1, 8, 22]),
            (0, [5, 9, 23]),
            (0, [7, 10, 19]),
            (0, [3, 15, 17]),
            (0, [2, 11, 18]),
            (0, [6, 10, 19]),
        ];

        for (center, pts) in arcs {
            utils::push_arc(vtx[center], vtx[pts[0]], vtx[pts[1]], nsubdivs, &mut vtx);
        }

        (vtx, idx)
    }
}
