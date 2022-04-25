use crate::math::Real;
use crate::shape::RoundCuboid;
use crate::transformation::utils;
use na::{self, Point2};

impl RoundCuboid {
    /// Discretize the boundary of this round cuboid as a polygonal line.
    pub fn to_polyline(&self, border_subdivs: u32) -> Vec<Point2<Real>> {
        let mut out_vtx = vec![];
        let he = self.inner_shape.half_extents;
        let br = self.border_radius;

        let arc_centers = [
            Point2::new(-he.x, -he.y),
            Point2::new(he.x, -he.y),
            Point2::new(he.x, he.y),
            Point2::new(-he.x, he.y),
        ];
        let arc_vertices = [
            (
                Point2::new(-he.x - br, -he.y),
                Point2::new(-he.x, -he.y - br),
            ),
            (Point2::new(he.x, -he.y - br), Point2::new(he.x + br, -he.y)),
            (Point2::new(he.x + br, he.y), Point2::new(he.x, he.y + br)),
            (Point2::new(-he.x, he.y + br), Point2::new(-he.x - br, he.y)),
        ];

        for i in 0..4 {
            out_vtx.push(arc_vertices[i].0);
            utils::push_arc(
                arc_centers[i],
                arc_vertices[i].0,
                arc_vertices[i].1,
                border_subdivs,
                &mut out_vtx,
            );
            out_vtx.push(arc_vertices[i].1);
        }

        out_vtx
    }
}
