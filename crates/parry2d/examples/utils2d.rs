use core::f32::consts::{FRAC_PI_2, FRAC_PI_4};
use std::sync::Arc;

use kiss3d::prelude::*;
use parry2d::shape::TriMesh;

/// As this file is used as a module from other examples,
/// rustc warns about dead code:
/// - `main()` is needed for this file to be included in examples
/// - For other functions, they may be "dead code" for an example, but not for others.
#[allow(dead_code)]
fn main() {
    println!(
        "This module contains helper functions to use kiss3d,
    isolated from the rest of the examples for the sake of simplicity."
    );
}

/// Uses [`kiss3d`] to display the line passed as parameter.
#[allow(dead_code)]
pub fn draw_polyline(window: &mut Window, polyline: Vec<(Vec2, Vec2)>, color: Color) {
    for line in polyline {
        let a = line.0;
        let b = line.1;
        draw_line_2d(window, a, b, color);
    }
}

/// Draws a text in the top left corner of the screen.
///
/// This uses a hardcoded position, size, color.
#[allow(dead_code)]
pub fn draw_text(window: &mut Window, font: &Arc<Font>, text: &str) {
    window.draw_text(text, Vec2::new(10.0, 66.0), 30.0, font, WHITE);
}

/// Returns [lissajous curve](https://en.wikipedia.org/wiki/Lissajous_curve) coordinates for time `t`.
///
/// This uses hardcoded parameters to have an arbitrary pleasing trajectory.
#[allow(dead_code)]
pub fn lissajous_2d(t: f32) -> Vec2 {
    // Some hardcoded parameters to have a pleasing lissajous trajectory.
    lissajous_2d_with_params(t, 3.0, 2.0, FRAC_PI_2, FRAC_PI_4)
}

/// Returns [lissajous curve](https://en.wikipedia.org/wiki/Lissajous_curve) coordinates.
#[allow(dead_code)]
pub fn lissajous_2d_with_params(t: f32, a: f32, b: f32, delta_x: f32, delta_y: f32) -> Vec2 {
    // Some hardcoded parameters to have a pleasing lissajous trajectory.

    let x = (a * t + delta_x).sin();
    let y = (b * t + delta_y).sin();
    Vec2::new(x, y) * 0.75f32
}

/// Uses [`kiss3d`] to display the line passed as parameter.
#[allow(dead_code)]
pub fn draw_line_2d(window: &mut Window, a: Vec2, b: Vec2, color: Color) {
    window.draw_line_2d(a, b, color, 2.0);
}

/// Uses [`kiss3d`] to display the trimesh passed as parameter.
#[allow(dead_code)]
pub fn draw_trimesh2(window: &mut Window, trimesh: &TriMesh, offset: Vec2) {
    let vertices = trimesh.vertices();
    for v in trimesh.indices() {
        let v0 = vertices[v[0] as usize] + offset;
        let v1 = vertices[v[1] as usize] + offset;
        let v2 = vertices[v[2] as usize] + offset;

        window.draw_line_2d(v0, v1, WHITE, 2.0);
        window.draw_line_2d(v0, v2, WHITE, 2.0);
        window.draw_line_2d(v2, v1, WHITE, 2.0);
    }
}

/// Uses [`kiss3d`] to display a wireframe of the polygon.
#[allow(dead_code)]
pub fn draw_polygon(window: &mut Window, polygon: &[Vec2], scale: f32, shift: Vec2, color: Color) {
    for i in 0..polygon.len() {
        let a = polygon[i];
        let b = polygon[(i + 1) % polygon.len()];
        window.draw_line_2d(a * scale + shift, b * scale + shift, color, 2.0);
    }
}

/// Uses [`kiss3d`] to display a cross, representing a point.
#[allow(dead_code)]
pub fn draw_point(window: &mut Window, point: Vec2, scale: f32, shift: Vec2, color: Color) {
    let edge_len = 0.15;
    let scaled_point = point * scale + shift;
    let edge_offset = Vec2::new(edge_len * scale, 0.0);
    let vert_offset = Vec2::new(0.0, edge_len * scale);
    window.draw_line_2d(
        scaled_point - edge_offset,
        scaled_point + edge_offset,
        color,
        2.0,
    );
    window.draw_line_2d(
        scaled_point - vert_offset,
        scaled_point + vert_offset,
        color,
        2.0,
    );
}

/// Draws a circle outline.
#[allow(dead_code)]
pub fn draw_circle(window: &mut Window, center: Vec2, radius: f32, color: Color) {
    let segments = 32;
    let tau = std::f32::consts::TAU;

    for i in 0..segments {
        let angle1 = (i as f32 / segments as f32) * tau;
        let angle2 = ((i + 1) as f32 / segments as f32) * tau;

        let p1 = center + Vec2::new(radius * angle1.cos(), radius * angle1.sin());
        let p2 = center + Vec2::new(radius * angle2.cos(), radius * angle2.sin());
        window.draw_line_2d(p1, p2, color, 2.0);
    }
}

/// Draws a 2D AABB (axis-aligned bounding box) as a rectangle outline.
#[allow(dead_code)]
pub fn draw_aabb2(window: &mut Window, mins: Vec2, maxs: Vec2, color: Color) {
    window.draw_line_2d(
        Vec2::new(mins.x, mins.y),
        Vec2::new(maxs.x, mins.y),
        color,
        2.0,
    );
    window.draw_line_2d(
        Vec2::new(maxs.x, mins.y),
        Vec2::new(maxs.x, maxs.y),
        color,
        2.0,
    );
    window.draw_line_2d(
        Vec2::new(maxs.x, maxs.y),
        Vec2::new(mins.x, maxs.y),
        color,
        2.0,
    );
    window.draw_line_2d(
        Vec2::new(mins.x, maxs.y),
        Vec2::new(mins.x, mins.y),
        color,
        2.0,
    );
}
