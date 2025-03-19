use crate::math::{Point, Vector};

/// Computes the index of the support point of a cloud of points.
#[inline]
pub fn point_cloud_support_point_id(dir: &Vector, points: &[Point]) -> usize {
    let mut best_pt = 0;
    let mut best_dot = points[0].coords.dot(dir);

    for (i, p) in points.iter().enumerate().skip(1) {
        let dot = p.coords.dot(dir);

        if dot > best_dot {
            best_dot = dot;
            best_pt = i;
        }
    }

    best_pt
}

/// Computes the support point of a cloud of points.
#[inline]
pub fn point_cloud_support_point(dir: &Vector, points: &[Point]) -> Point {
    points[point_cloud_support_point_id(dir, points)]
}
