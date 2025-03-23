//! Function to check if a point is inside a triangle and related functions.

use crate::math::{Point, Real};

#[derive(Eq, PartialEq, Debug, Copy, Clone)]
/// The orientation or winding direction of a corner or polygon.
pub enum Orientation {
    /// Counter-clockwise
    Ccw,
    /// Clockwise
    Cw,
    /// Neither (a straight line)
    None,
}

/// Returns the direction of a line through `p1`, `p2` and `p3`.
///
/// Counter-clockwise example:
/// o p1
///  .        o p3
///   .     .
///    .  .
///     o p2
///
/// Clockwise example:
///     o p2
///    .  .
///   .     .
///  .        o p3
/// o p1
pub fn corner_direction(p1: &Point<Real>, p2: &Point<Real>, p3: &Point<Real>) -> Orientation {
    let v1 = p1 - p2;
    let v2 = p3 - p2;
    let cross: Real = v1.perp(&v2);

    match cross
        .partial_cmp(&0.0)
        .expect("Found NaN while computing corner direction.")
    {
        core::cmp::Ordering::Less => Orientation::Ccw,
        core::cmp::Ordering::Equal => Orientation::None,
        core::cmp::Ordering::Greater => Orientation::Cw,
    }
}

/// Returns `true` if point `p` is in triangle with corners `v1`, `v2` and `v3`.
/// Returns `None` if the triangle is invalid i.e. all points are the same or on a straight line.
pub fn is_point_in_triangle(
    p: &Point<Real>,
    v1: &Point<Real>,
    v2: &Point<Real>,
    v3: &Point<Real>,
) -> Option<bool> {
    let d1 = corner_direction(p, v1, v2);
    let d2 = corner_direction(p, v2, v3);
    let d3 = corner_direction(p, v3, v1);

    let has_cw = d1 == Orientation::Cw || d2 == Orientation::Cw || d3 == Orientation::Cw;
    let has_ccw = d1 == Orientation::Ccw || d2 == Orientation::Ccw || d3 == Orientation::Ccw;

    if d1 == Orientation::None && d2 == Orientation::None && d3 == Orientation::None {
        None
    } else {
        Some(!(has_cw && has_ccw))
    }
}
