use crate::bounding_volume;
use crate::bounding_volume::BoundingSphere;
use crate::math::{Isometry, Real};
use crate::shape::ConvexPolyhedron;

impl ConvexPolyhedron {
    #[inline]
    pub fn bounding_sphere(&self, m: &Isometry<Real>) -> BoundingSphere {
        let bv: BoundingSphere = self.local_bounding_sphere();
        bv.transform_by(m)
    }

    #[inline]
    pub fn local_bounding_sphere(&self) -> BoundingSphere {
        let (center, radius) = bounding_volume::details::point_cloud_bounding_sphere(self.points());

        BoundingSphere::new(center, radius)
    }
}
