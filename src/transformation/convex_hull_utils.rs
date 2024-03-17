use crate::math::*;
use crate::num::Bounded;
#[cfg(feature = "dim3")]
use {crate::bounding_volume, crate::math::Point};

/// Returns the index of the support point of a list of points.
pub fn support_point_id(direction: &Vector, points: &[Point]) -> Option<usize> {
    let mut argmax = None;
    let _max: Real = Bounded::max_value();
    let mut max = -_max;

    for (id, pt) in points.iter().enumerate() {
        let dot = direction.dot(pt.as_vector());

        if dot > max {
            argmax = Some(id);
            max = dot;
        }
    }

    argmax
}

/// Returns the index of the support point of a list of 2D points.
pub fn support_point_id_2d(direction: &Vector2, points: &[Point2]) -> Option<usize> {
    let mut argmax = None;
    let _max: Real = Bounded::max_value();
    let mut max = -_max;

    for (id, pt) in points.iter().enumerate() {
        let dot = direction.dot(pt.as_vector());

        if dot > max {
            argmax = Some(id);
            max = dot;
        }
    }

    argmax
}

/// Returns the index of the support point of an indexed list of points.
pub fn indexed_support_point_id<I>(direction: &Vector, points: &[Point], idx: I) -> Option<usize>
where
    I: Iterator<Item = usize>,
{
    let mut argmax = None;
    let mut max = -Real::MAX;

    for i in idx.into_iter() {
        let dot = direction.dot(points[i].as_vector());

        if dot > max {
            argmax = Some(i);
            max = dot;
        }
    }

    argmax
}

/// Returns the index of the support point of an indexed list of 2D points.
pub fn indexed_support_point_id_2d<I>(
    direction: &Vector2,
    points: &[Point2],
    idx: I,
) -> Option<usize>
where
    I: Iterator<Item = usize>,
{
    let mut argmax = None;
    let mut max = -Real::MAX;

    for i in idx.into_iter() {
        let dot = direction.dot(points[i].as_vector());

        if dot > max {
            argmax = Some(i);
            max = dot;
        }
    }

    argmax
}

/// Returns the number `n` such that `points[idx.nth(n)]` is the support point.
#[cfg(feature = "dim3")] // We only use this in 3D right now.
pub fn indexed_support_point_nth<I>(direction: &Vector, points: &[Point], idx: I) -> Option<usize>
where
    I: Iterator<Item = usize>,
{
    let mut argmax = None;
    let mut max = -Real::MAX;

    for (k, i) in idx.into_iter().enumerate() {
        let dot = direction.dot(points[i].as_vector());

        if dot > max {
            argmax = Some(k);
            max = dot;
        }
    }

    argmax
}

/// Scale and center the given set of point depending on their Aabb.
#[cfg(feature = "dim3")]
pub fn normalize(coords: &mut [Point]) -> (Point, Real) {
    let aabb = bounding_volume::details::local_point_cloud_aabb(&coords[..]);
    let diag = distance(aabb.mins, aabb.maxs);
    let center = aabb.center();

    for c in coords.iter_mut() {
        *c = (*c + (-center.as_vector())) / diag;
    }

    (center, diag)
}
