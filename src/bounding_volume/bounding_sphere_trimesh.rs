use crate::bounding_volume::BoundingSphere;
use crate::math::{Isometry, Real};
use crate::shape::TriMesh;

impl TriMesh {
    #[inline]
    pub fn bounding_sphere(&self, m: &Isometry<Real>) -> BoundingSphere {
        self.local_aabb().bounding_sphere().transform_by(m)
    }

    #[inline]
    pub fn local_bounding_sphere(&self) -> BoundingSphere {
        self.local_aabb().bounding_sphere()
    }
}
