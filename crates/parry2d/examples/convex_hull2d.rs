mod common_macroquad2d;

use core::f32::consts::{FRAC_PI_2, FRAC_PI_4};

use common_macroquad2d::{draw_point, draw_polygon, lissajous_2d_with_params, na_from_mquad};
use macroquad::prelude::*;
use parry2d::math::Vector;
use parry2d::transformation;

const RENDER_SCALE: f32 = 30.0;

#[macroquad::main("convex_hull2d")]
async fn main() {
    let count = 9;
    let mut pts = vec![Vector::ZERO; count];

    let render_pos = Vector::new(300.0, 300.0);

    loop {
        let elapsed_time = get_time() as f32;
        let elapsed_time_slow = elapsed_time * 0.2;
        clear_background(BLACK);

        for (i, pt) in pts.iter_mut().enumerate() {
            *pt = na_from_mquad(lissajous_2d_with_params(
                (i * i) as f32 + elapsed_time_slow,
                2.0 + i as f32 / 3.0,
                (i as f32 / count as f32) + elapsed_time_slow.cos() * 0.1,
                (elapsed_time_slow as f32 + i as f32).cos() * 0.1 + FRAC_PI_2,
                FRAC_PI_4,
            )) * 5f32;
            draw_point(*pt, RENDER_SCALE, render_pos, RED);
        }

        /*
         *
         * Compute the convex hull.
         *
         */
        let convex_hull = transformation::convex_hull(&pts);
        draw_polygon(&convex_hull, RENDER_SCALE, render_pos, WHITE);
        next_frame().await
    }
}
