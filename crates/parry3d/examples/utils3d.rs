use core::f32::consts::{FRAC_PI_2, FRAC_PI_4, FRAC_PI_6};
use std::sync::Arc;

use kiss3d::prelude::*;
use kiss3d::procedural::{IndexBuffer, RenderMesh};

#[allow(dead_code)]
fn main() {
    println!(
        "This module contains helper functions to use kiss3d,
    isolated from the rest of the examples for the sake of simplicity."
    );
}

/// Converts a hue (from 0..=1) to rgb
#[allow(dead_code)]
pub fn hue_to_rgb(h: f32) -> (f32, f32, f32) {
    let kr = (5.0 + h * 6.0).rem_euclid(6.0);
    let kg = (3.0 + h * 6.0).rem_euclid(6.0);
    let kb = (1.0 + h * 6.0).rem_euclid(6.0);

    let r = 1.0 - kr.min(4.0 - kr).clamp(0.0, 1.0);
    let g = 1.0 - kg.min(4.0 - kg).clamp(0.0, 1.0);
    let b = 1.0 - kb.min(4.0 - kb).clamp(0.0, 1.0);

    (r, g, b)
}

/// Returns [lissajous curve](https://en.wikipedia.org/wiki/Lissajous_curve) coordinates for time `t`.
///
/// This uses hardcoded parameters to have an arbitrary pleasing trajectory.
#[allow(dead_code)]
pub fn lissajous_3d(t: f32) -> Vec3 {
    // Some hardcoded parameters to have a pleasing lissajous trajectory.
    lissajous_3d_with_params(t, 3.0, 2.0, 1.0, FRAC_PI_2, FRAC_PI_4, FRAC_PI_6)
}

/// Returns [lissajous curve](https://en.wikipedia.org/wiki/Lissajous_curve) coordinates.
#[allow(dead_code)]
pub fn lissajous_3d_with_params(
    t: f32,
    a: f32,
    b: f32,
    c: f32,
    delta_x: f32,
    delta_y: f32,
    delta_z: f32,
) -> Vec3 {
    let x = (a * t + delta_x).sin();
    let y = (b * t + delta_y).sin();
    let z = (c * t + delta_z).sin();
    Vec3::new(x, y, z) * 0.75f32
}

/// Uses [`kiss3d`] to display the line passed as parameter.
#[allow(dead_code)]
pub fn draw_polyline(window: &mut Window, polyline: Vec<(Vec3, Vec3)>, color: Color) {
    for line in polyline {
        let a = line.0;
        let b = line.1;
        window.draw_line(a, b, color, 2.0, false);
    }
}

/// Draws a text in the top left corner of the screen.
///
/// This uses a hardcoded position, size, color.
#[allow(dead_code)]
pub fn draw_text(window: &mut Window, font: &Arc<Font>, text: &str) {
    window.draw_text(text, Vec2::new(10.0, 66.0), 30.0, font, WHITE);
}

/// Creates a GpuMesh3d from parry3d's trimesh vertices and indices.
#[allow(dead_code)]
pub fn create_mesh_from_trimesh(
    vertices: Vec<Vec3>,
    indices: Vec<[u32; 3]>,
) -> std::rc::Rc<std::cell::RefCell<kiss3d::resource::GpuMesh3d>> {
    use std::cell::RefCell;
    use std::rc::Rc;

    let mut render_mesh =
        RenderMesh::new(vertices, None, None, Some(IndexBuffer::Unified(indices)));
    render_mesh.replicate_vertices();
    render_mesh.recompute_normals();
    Rc::new(RefCell::new(render_mesh.into()))
}

/// Draws a 3D AABB (axis-aligned bounding box) as a wireframe.
#[allow(dead_code)]
pub fn draw_aabb3(window: &mut Window, mins: Vec3, maxs: Vec3, color: Color) {
    // Bottom face
    window.draw_line(
        Vec3::new(mins.x, mins.y, mins.z),
        Vec3::new(maxs.x, mins.y, mins.z),
        color,
        2.0,
        false,
    );
    window.draw_line(
        Vec3::new(maxs.x, mins.y, mins.z),
        Vec3::new(maxs.x, mins.y, maxs.z),
        color,
        2.0,
        false,
    );
    window.draw_line(
        Vec3::new(maxs.x, mins.y, maxs.z),
        Vec3::new(mins.x, mins.y, maxs.z),
        color,
        2.0,
        false,
    );
    window.draw_line(
        Vec3::new(mins.x, mins.y, maxs.z),
        Vec3::new(mins.x, mins.y, mins.z),
        color,
        2.0,
        false,
    );
    // Top face
    window.draw_line(
        Vec3::new(mins.x, maxs.y, mins.z),
        Vec3::new(maxs.x, maxs.y, mins.z),
        color,
        2.0,
        false,
    );
    window.draw_line(
        Vec3::new(maxs.x, maxs.y, mins.z),
        Vec3::new(maxs.x, maxs.y, maxs.z),
        color,
        2.0,
        false,
    );
    window.draw_line(
        Vec3::new(maxs.x, maxs.y, maxs.z),
        Vec3::new(mins.x, maxs.y, maxs.z),
        color,
        2.0,
        false,
    );
    window.draw_line(
        Vec3::new(mins.x, maxs.y, maxs.z),
        Vec3::new(mins.x, maxs.y, mins.z),
        color,
        2.0,
        false,
    );
    // Vertical edges
    window.draw_line(
        Vec3::new(mins.x, mins.y, mins.z),
        Vec3::new(mins.x, maxs.y, mins.z),
        color,
        2.0,
        false,
    );
    window.draw_line(
        Vec3::new(maxs.x, mins.y, mins.z),
        Vec3::new(maxs.x, maxs.y, mins.z),
        color,
        2.0,
        false,
    );
    window.draw_line(
        Vec3::new(maxs.x, mins.y, maxs.z),
        Vec3::new(maxs.x, maxs.y, maxs.z),
        color,
        2.0,
        false,
    );
    window.draw_line(
        Vec3::new(mins.x, mins.y, maxs.z),
        Vec3::new(mins.x, maxs.y, maxs.z),
        color,
        2.0,
        false,
    );
}

/// Draws a sphere wireframe using 3 orthogonal circles.
#[allow(dead_code)]
pub fn draw_sphere_wires(window: &mut Window, center: Vec3, radius: f32, color: Color) {
    let segments = 32;
    let tau = std::f32::consts::TAU;

    for i in 0..segments {
        let angle1 = (i as f32 / segments as f32) * tau;
        let angle2 = ((i + 1) as f32 / segments as f32) * tau;

        // XY plane circle
        let p1 = center + Vec3::new(radius * angle1.cos(), radius * angle1.sin(), 0.0);
        let p2 = center + Vec3::new(radius * angle2.cos(), radius * angle2.sin(), 0.0);
        window.draw_line(p1, p2, color, 2.0, false);

        // XZ plane circle
        let p1 = center + Vec3::new(radius * angle1.cos(), 0.0, radius * angle1.sin());
        let p2 = center + Vec3::new(radius * angle2.cos(), 0.0, radius * angle2.sin());
        window.draw_line(p1, p2, color, 2.0, false);

        // YZ plane circle
        let p1 = center + Vec3::new(0.0, radius * angle1.cos(), radius * angle1.sin());
        let p2 = center + Vec3::new(0.0, radius * angle2.cos(), radius * angle2.sin());
        window.draw_line(p1, p2, color, 2.0, false);
    }
}
