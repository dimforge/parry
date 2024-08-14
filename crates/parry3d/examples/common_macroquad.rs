use std::f32::consts::{FRAC_PI_2, FRAC_PI_4, FRAC_PI_6};

use macroquad::{
    color::{Color, WHITE},
    math::Vec3,
    models::draw_line_3d,
};
use nalgebra::Point3;
use parry3d::math::Real;

fn main() {
    println!(
        "This module contains helper fubnctions to use macroquad,
    isolated from the rest of the examples for the sake of simplicity."
    );
}

pub fn mquad_from_na(a: Point3<Real>) -> Vec3 {
    Vec3::new(a.x, a.y, a.z)
}

pub fn na_from_mquad(a: Vec3) -> Point3<Real> {
    Point3::new(a.x, a.y, a.z)
}

pub fn lissajous_3d(t: f32) -> Vec3 {
    // Some hardcoded parameters to have a pleasing lissajous trajectory.
    let (a, b, c, delta_x, delta_y, delta_z) = (3.0, 2.0, 1.0, FRAC_PI_2, FRAC_PI_4, FRAC_PI_6);

    let x = (a * t + delta_x).sin();
    let y = (b * t + delta_y).sin();
    let z = (c * t + delta_z).sin();
    Vec3::new(x, y, z) * 0.75f32
}

pub fn draw_polyline(polygon: Vec<(Vec3, Vec3)>, color: Color) {
    for i in 0..polygon.len() {
        let a = polygon[i].0;
        let b = polygon[i].1;
        draw_line_3d(a, b, color);
    }
}

pub fn easy_draw_text(text: &str) {
    macroquad::text::draw_text(text, 10.0, 48.0 + 18.0, 30.0, WHITE);
}
