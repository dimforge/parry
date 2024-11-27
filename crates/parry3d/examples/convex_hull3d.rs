mod common_macroquad3d;

use std::f32::consts::{FRAC_PI_2, FRAC_PI_4, FRAC_PI_6};

use common_macroquad3d::{
    lissajous_3d_with_params, mquad_from_na, mquad_mesh_from_points, na_from_mquad,
};
use macroquad::prelude::*;
use nalgebra::Point3;
use parry3d::transformation;

#[macroquad::main("convex_hull3d")]
async fn main() {
    let count = 9;
    let mut pts = vec![Point3::default(); count];

    let camera_pos = Vec3::new(8.0, 8.0, 8.0);
    loop {
        let elapsed_time = get_time() as f32;
        let elapsed_time_slow = elapsed_time * 0.2;
        clear_background(BLACK);

        for (i, pt) in pts.iter_mut().enumerate() {
            *pt = na_from_mquad(lissajous_3d_with_params(
                (i * i) as f32 + elapsed_time_slow,
                2.0 + i as f32 / 3.0,
                1f32 + (i as f32).sin() * 0.2,
                (i as f32 / count as f32) + elapsed_time_slow.cos() * 0.1,
                (elapsed_time_slow as f32 + i as f32).cos() * 0.1 + FRAC_PI_2,
                FRAC_PI_4,
                FRAC_PI_6,
            )) * 5f32;
            draw_sphere(mquad_from_na(*pt), 0.1f32, None, RED);
        }
        // Initialize 3D camera.
        set_camera(&Camera3D {
            position: camera_pos,
            up: Vec3::new(0f32, 1f32, 0f32),
            target: Vec3::new(0.5f32, 0f32, 0.5f32),
            ..Default::default()
        });
        /*
         *
         * Compute the convex hull.
         *
         */
        let convex_hull = transformation::convex_hull(&pts);
        let mesh = mquad_mesh_from_points(&convex_hull, Vec3::new(5.0, 10.0, 3.0), DARKGRAY);
        draw_mesh(&mesh);
        next_frame().await
    }
}
