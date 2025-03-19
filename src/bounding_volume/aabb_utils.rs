use std::iter::IntoIterator;

use crate::bounding_volume::Aabb;
use crate::math::{Isometry, Point, Vector, DIM};
use crate::shape::SupportMap;
use na;

/// Computes the [`Aabb`] of an [support mapped shape](SupportMap).
#[cfg(feature = "dim3")]
pub fn support_map_aabb<G>(m: &Isometry, i: &G) -> Aabb
where
    G: SupportMap,
{
    let mut min = na::zero::<Vector>();
    let mut max = na::zero::<Vector>();
    let mut basis = na::zero::<Vector>();

    for d in 0..DIM {
        // TODO: this could be further improved iterating on `m`'s columns, and passing
        // Id as the transformation matrix.
        basis[d] = 1.0;
        max[d] = i.support_point(m, &basis)[d];

        basis[d] = -1.0;
        min[d] = i.support_point(m, &basis)[d];

        basis[d] = 0.0;
    }

    Aabb::new(Point::from(min), Point::from(max))
}

/// Computes the [`Aabb`] of an [support mapped shape](SupportMap).
pub fn local_support_map_aabb<G>(i: &G) -> Aabb
where
    G: SupportMap,
{
    let mut min = na::zero::<Vector>();
    let mut max = na::zero::<Vector>();
    let mut basis = na::zero::<Vector>();

    for d in 0..DIM {
        // TODO: this could be further improved iterating on `m`'s columns, and passing
        // Id as the transformation matrix.
        basis[d] = 1.0;
        max[d] = i.local_support_point(&basis)[d];

        basis[d] = -1.0;
        min[d] = i.local_support_point(&basis)[d];

        basis[d] = 0.0;
    }

    Aabb::new(Point::from(min), Point::from(max))
}

/// Computes the [`Aabb`] of a set of points transformed by `m`.
pub fn point_cloud_aabb<'a, I>(m: &Isometry, pts: I) -> Aabb
where
    I: IntoIterator<Item = &'a Point>,
{
    let mut it = pts.into_iter();

    let p0 = it.next().expect(
        "Point cloud Aabb construction: the input iterator should yield at least one point.",
    );
    let wp0 = m.transform_point(p0);
    let mut min: Point = wp0;
    let mut max: Point = wp0;

    for pt in it {
        let wpt = m * pt;
        min = min.inf(&wpt);
        max = max.sup(&wpt);
    }

    Aabb::new(min, max)
}

/// Computes the [`Aabb`] of a set of points.
pub fn local_point_cloud_aabb<'a, I>(pts: I) -> Aabb
where
    I: IntoIterator<Item = &'a Point>,
{
    let mut it = pts.into_iter();

    let p0 = it.next().expect(
        "Point cloud Aabb construction: the input iterator should yield at least one point.",
    );
    let mut min: Point = *p0;
    let mut max: Point = *p0;

    for pt in it {
        min = min.inf(pt);
        max = max.sup(pt);
    }

    Aabb::new(min, max)
}
