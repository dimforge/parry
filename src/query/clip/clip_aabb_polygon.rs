use crate::bounding_volume::AABB;
use crate::math::{Point, Real, Vector};

impl AABB {
    /// Computes the intersections between this AABB and the given polygon.
    ///
    /// The results is written into `points` directly. The input points are
    /// assumed to form a convex polygon where all points lie on the same plane.
    /// In order to avoid internal allocations, uses `self.clip_polygon_with_workspace`
    /// instead.
    #[inline]
    pub fn clip_polygon(&self, points: &mut Vec<Point<Real>>) {
        let mut workspace = Vec::new();
        self.clip_polygon_with_workspace(points, &mut workspace)
    }

    /// Computes the intersections between this AABB and the given polygon.
    ///
    /// The results is written into `points` directly. The input points are
    /// assumed to form a convex polygon where all points lie on the same plane.
    #[inline]
    pub fn clip_polygon_with_workspace(
        &self,
        points: &mut Vec<Point<Real>>,
        workspace: &mut Vec<Point<Real>>,
    ) {
        super::clip_halfspace_polygon(&self.mins, &-Vector::x(), &points, workspace);
        super::clip_halfspace_polygon(&self.maxs, &Vector::x(), &workspace, points);

        super::clip_halfspace_polygon(&self.mins, &-Vector::y(), &points, workspace);
        super::clip_halfspace_polygon(&self.maxs, &Vector::y(), &workspace, points);

        #[cfg(feature = "dim3")]
        {
            super::clip_halfspace_polygon(&self.mins, &-Vector::z(), &points, workspace);
            super::clip_halfspace_polygon(&self.maxs, &Vector::z(), &workspace, points);
        }
    }
}
