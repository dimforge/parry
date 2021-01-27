use crate::math::Real;
use crate::shape::HeightField;
use na::Point3;

impl HeightField {
    /// Discretize the boundary of this ball as a triangle-mesh.
    pub fn to_trimesh(&self) -> (Vec<Point3<Real>>, Vec<[u32; 3]>) {
        let mut vertices = Vec::new();
        let mut indices = Vec::new();

        for (i, tri) in self.triangles().enumerate() {
            vertices.push(tri.a);
            vertices.push(tri.b);
            vertices.push(tri.c);

            let i = i as u32;
            indices.push([i * 3, i * 3 + 1, i * 3 + 2])
        }

        (vertices, indices)
    }
}
