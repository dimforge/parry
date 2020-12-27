use crate::shape::ConvexPolyhedron;
use na::Point3;

impl ConvexPolyhedron {
    pub fn to_trimesh(&self) -> (Vec<Point3<f32>>, Vec<Point3<u32>>) {
        let mut indices = Vec::new();

        for face in self.faces() {
            let i1 = face.first_vertex_or_edge;
            let i2 = i1 + face.num_vertices_or_edges;
            let first_id = self.vertices_adj_to_face()[i1] as u32;

            for idx in self.vertices_adj_to_face()[i1 + 1..i2].windows(2) {
                indices.push(Point3::new(first_id, idx[0] as u32, idx[1] as u32));
            }
        }

        (self.points().to_vec(), indices)
    }
}
