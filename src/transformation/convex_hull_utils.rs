use crate::bounding_volume;
use crate::math::{Point, Real};
use crate::num::Bounded;
use na;
use na::allocator::Allocator;
use na::base::{DefaultAllocator, DimName};

/// Returns the index of the support point of a list of points.
pub fn support_point_id<D: DimName>(
    direction: &na::VectorN<Real, D>,
    points: &[na::Point<Real, D>],
) -> Option<usize>
where
    DefaultAllocator: Allocator<Real, D>,
{
    let mut argmax = None;
    let _max: Real = Bounded::max_value();
    let mut max = -_max;

    for (id, pt) in points.iter().enumerate() {
        let dot = direction.dot(&pt.coords);

        if dot > max {
            argmax = Some(id);
            max = dot;
        }
    }

    argmax
}

/// Returns the index of the support point of an indexed list of points.
pub fn indexed_support_point_id<D: DimName>(
    direction: &na::VectorN<Real, D>,
    points: &[na::Point<Real, D>],
    idx: &[usize],
) -> Option<usize>
where
    DefaultAllocator: Allocator<Real, D>,
{
    let mut argmax = None;
    let _max: Real = Bounded::max_value();
    let mut max = -_max;

    for i in idx.iter() {
        let dot = direction.dot(&points[*i].coords);

        if dot > max {
            argmax = Some(*i);
            max = dot;
        }
    }

    argmax
}

/// Scale and center the given set of point depending on their AABB.
pub fn normalize(coords: &mut [Point<Real>]) -> (Point<Real>, Real) {
    let aabb = bounding_volume::local_point_cloud_aabb(&coords[..]);
    let diag = na::distance(&aabb.mins, &aabb.maxs);
    let center = aabb.center();

    for c in coords.iter_mut() {
        *c = (*c + (-center.coords)) / diag;
    }

    (center, diag)
}

/// Scale and translates the given set of point.
pub fn denormalize(coords: &mut [Point<Real>], center: &Point<Real>, diag: Real) {
    for c in coords.iter_mut() {
        *c = *c * diag + center.coords;
    }
}
