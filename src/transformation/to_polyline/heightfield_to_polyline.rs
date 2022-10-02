use crate::math::Real;
use crate::shape::HeightField;
use na::Point2;

impl HeightField {
    /// Rasterize this heightfield as a (potentially discontinuous) polyline.
    pub fn to_polyline(&self) -> (Vec<Point2<Real>>, Vec<[u32; 2]>) {
        let mut vertices = vec![];
        let mut indices = vec![];

        for seg in self.segments() {
            let base_id = vertices.len() as u32;
            if vertices.last().map(|pt| seg.a != *pt).unwrap_or(true) {
                indices.push([base_id, base_id + 1]);
                vertices.push(seg.a);
            } else {
                indices.push([base_id - 1, base_id]);
            }

            vertices.push(seg.b);
        }

        (vertices, indices)
    }
}
