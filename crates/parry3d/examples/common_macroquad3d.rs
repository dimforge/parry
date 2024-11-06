use std::f32::consts::{FRAC_PI_2, FRAC_PI_4, FRAC_PI_6};

use macroquad::{
    color::{Color, WHITE},
    math::{Vec2, Vec3, Vec4},
    models::{draw_line_3d, Mesh},
    ui::Vertex,
};
use nalgebra::Point3;
use parry3d::math::Real;

#[allow(dead_code)]
fn main() {
    println!(
        "This module contains helper functions to use macroquad,
    isolated from the rest of the examples for the sake of simplicity."
    );
}

/// Converts a [`nalgebra::Point3`] to a [`Vec3`], which is used by [`macroquad`]
#[allow(dead_code)]
pub fn mquad_from_na(a: Point3<Real>) -> Vec3 {
    Vec3::new(a.x, a.y, a.z)
}

/// Converts a [`Vec3`] to a [`nalgebra::Point3`], which is used by [`parry3d`]
#[allow(dead_code)]
pub fn na_from_mquad(a: Vec3) -> Point3<Real> {
    Point3::new(a.x, a.y, a.z)
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

/// Uses [`macroquad`] to display the line passed as parameter.
#[allow(dead_code)]
pub fn draw_polyline(polyline: Vec<(Vec3, Vec3)>, color: Color) {
    for line in polyline {
        let a = line.0;
        let b = line.1;
        draw_line_3d(a, b, color);
    }
}

/// Draws a text in the top left corner of the screen.
///
/// This uses a hardcoded position, size, color.
#[allow(dead_code)]
pub fn easy_draw_text(text: &str) {
    macroquad::text::draw_text(text, 10.0, 48.0 + 18.0, 30.0, WHITE);
}

/// Create a usable mesh for [`macroquad`].
///
/// This duplicates the trimesh vertices, computes their normals,
/// and bakes light into its vertices colors using [`mquad_compute_normals_and_bake_light`].
#[allow(dead_code)]
pub fn mquad_mesh_from_points(
    trimesh: &(Vec<Point3<Real>>, Vec<[u32; 3]>),
    light_pos: Vec3,
    color: Color,
) -> Mesh {
    let (points, indices) = trimesh;
    // Transform the parry mesh into a mquad Mesh
    let (mquad_points, mquad_indices) = (
        points
            .iter()
            .map(|p| Vertex {
                position: mquad_from_na(*p),
                uv: Vec2::new(p.x, p.y),
                color: color.into(),
                normal: Vec4::ZERO,
            })
            .collect(),
        indices.iter().flatten().map(|v| *v as u16).collect(),
    );

    // Macroquad does support adding normals to vertices, but we´d have to provide shaders for them.
    // so we're baking a color into these vertices.
    // See https://github.com/not-fl3/macroquad/issues/321.

    // Compute the normal of each vertex, making them unique
    let vertices: Vec<Vertex> =
        mquad_compute_normals_and_bake_light(&mquad_points, &mquad_indices, light_pos);
    // Regenerate the index for each vertex.
    let indices: Vec<u16> = (0..vertices.len() * 3)
        .into_iter()
        .map(|i| i as u16)
        .collect();
    let mesh = Mesh {
        vertices,
        indices,
        texture: None,
    };
    mesh
}

/// Bakes light into vertices, using an hardcoded light strength.
#[allow(dead_code)]
pub fn mquad_compute_normals_and_bake_light(
    points: &Vec<Vertex>,
    indices: &Vec<u16>,
    light_pos: Vec3,
) -> Vec<Vertex> {
    let mut vertices: Vec<Vertex> = Vec::<Vertex>::new();
    for indices in indices.chunks(3) {
        let v0 = &points[indices[0] as usize];
        let v1 = &points[indices[1] as usize];
        let v2 = &points[indices[2] as usize];
        let normal = (v0.position - v2.position)
            .cross(v1.position - v2.position)
            .normalize();
        let brightness_mod = 0.4 + (0.6 / 2.) * (normal.dot(light_pos) + 1.);

        for &i in indices.iter() {
            let mut color = points[i as usize].color;
            color[0] = (color[0] as f32 * brightness_mod) as u8;
            color[1] = (color[1] as f32 * brightness_mod) as u8;
            color[2] = (color[2] as f32 * brightness_mod) as u8;

            vertices.push(Vertex {
                position: points[i as usize].position,
                uv: Vec2::ZERO,
                color: color,
                normal: Vec4::ZERO,
            });
        }
    }
    vertices
}
