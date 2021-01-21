//! Utilities useful for various generations tasks.

use crate::math::{Isometry, Point, Real, Vector};
use crate::na::ComplexField;
#[cfg(feature = "dim3")]
use {crate::math::DIM, num::Zero};

pub fn transformed(mut points: Vec<Point<Real>>, m: Isometry<Real>) -> Vec<Point<Real>> {
    points.iter_mut().for_each(|p| *p = m * *p);
    points
}

pub fn scaled(mut points: Vec<Point<Real>>, scale: Vector<Real>) -> Vec<Point<Real>> {
    points
        .iter_mut()
        .for_each(|p| p.coords.component_mul_assign(&scale));
    points
}

// FIXME: remove that in favor of `push_xy_circle` ?
/// Pushes a discretized counterclockwise circle to a buffer.
#[cfg(feature = "dim3")]
#[inline]
pub fn push_circle(radius: Real, nsubdiv: u32, dtheta: Real, y: Real, out: &mut Vec<Point<Real>>) {
    let mut curr_theta = Real::zero();

    for _ in 0..nsubdiv {
        out.push(Point::new(
            ComplexField::cos(curr_theta) * radius,
            y.clone(),
            ComplexField::sin(curr_theta) * radius,
        ));
        curr_theta = curr_theta + dtheta;
    }
}

/// Pushes a discretized counterclockwise circle to a buffer.
/// The circle is contained on the plane spanned by the `x` and `y` axis.
#[inline]
#[cfg(feature = "dim2")]
pub fn push_xy_arc(radius: Real, nsubdiv: u32, dtheta: Real, out: &mut Vec<Point<Real>>) {
    let mut curr_theta: Real = 0.0;

    for _ in 0..nsubdiv {
        let mut pt_coords = Vector::zeros();

        pt_coords[0] = ComplexField::cos(curr_theta) * radius;
        pt_coords[1] = ComplexField::sin(curr_theta) * radius;
        out.push(Point::from(pt_coords));

        curr_theta = curr_theta + dtheta;
    }
}

/// Creates the faces from two circles with the same discretization.
#[cfg(feature = "dim3")]
#[inline]
pub fn push_ring_indices(
    base_lower_circle: u32,
    base_upper_circle: u32,
    nsubdiv: u32,
    out: &mut Vec<[u32; DIM]>,
) {
    push_open_ring_indices(base_lower_circle, base_upper_circle, nsubdiv, out);

    // adjust the last two triangles
    push_rectangle_indices(
        base_upper_circle,
        base_upper_circle + nsubdiv - 1,
        base_lower_circle,
        base_lower_circle + nsubdiv - 1,
        out,
    );
}

/// Creates the faces from two circles with the same discretization.
#[cfg(feature = "dim3")]
#[inline]
pub fn push_open_ring_indices(
    base_lower_circle: u32,
    base_upper_circle: u32,
    nsubdiv: u32,
    out: &mut Vec<[u32; DIM]>,
) {
    assert!(nsubdiv > 0);

    for i in 0..nsubdiv - 1 {
        let bli = base_lower_circle + i;
        let bui = base_upper_circle + i;
        push_rectangle_indices(bui + 1, bui, bli + 1, bli, out);
    }
}

/// Creates the faces from a circle and a point that is shared by all triangle.
#[cfg(feature = "dim3")]
#[inline]
pub fn push_degenerate_top_ring_indices(
    base_circle: u32,
    point: u32,
    nsubdiv: u32,
    out: &mut Vec<[u32; DIM]>,
) {
    push_degenerate_open_top_ring_indices(base_circle, point, nsubdiv, out);

    out.push([base_circle + nsubdiv - 1, point, base_circle]);
}

/// Creates the faces from a circle and a point that is shared by all triangle.
#[cfg(feature = "dim3")]
#[inline]
pub fn push_degenerate_open_top_ring_indices(
    base_circle: u32,
    point: u32,
    nsubdiv: u32,
    out: &mut Vec<[u32; DIM]>,
) {
    assert!(nsubdiv > 0);

    for i in 0..nsubdiv - 1 {
        out.push([base_circle + i, point, base_circle + i + 1]);
    }
}

/// Pushes indices so that a circle is filled with triangles. Each triangle will have the
/// `base_circle` point in common.
/// Pushes `nsubdiv - 2` elements to `out`.
#[cfg(feature = "dim3")]
#[inline]
pub fn push_filled_circle_indices(base_circle: u32, nsubdiv: u32, out: &mut Vec<[u32; DIM]>) {
    for i in base_circle + 1..base_circle + nsubdiv - 1 {
        out.push([base_circle, i, i + 1]);
    }
}

/// Given four corner points, pushes to two counterclockwise triangles to `out`.
///
/// # Arguments:
/// * `ul` - the up-left point.
/// * `dl` - the down-left point.
/// * `dr` - the down-left point.
/// * `ur` - the up-left point.
#[cfg(feature = "dim3")]
#[inline]
pub fn push_rectangle_indices(ul: u32, ur: u32, dl: u32, dr: u32, out: &mut Vec<[u32; DIM]>) {
    out.push([ul.clone(), dl, dr.clone()]);
    out.push([dr, ur, ul]);
}

/// Reverses the clockwising of a set of faces.
#[cfg(feature = "dim3")]
#[inline]
pub fn reverse_clockwising(indices: &mut [[u32; DIM]]) {
    for i in indices.iter_mut() {
        i.swap(0, 1);
    }
}
