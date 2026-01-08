mod utils2d;

use core::f32::consts::{FRAC_PI_2, FRAC_PI_4};

use kiss3d::prelude::*;
use parry2d::transformation;
use utils2d::{draw_point, draw_polygon, lissajous_2d_with_params};

const RENDER_SCALE: f32 = 30.0;

#[kiss3d::main]
async fn main() {
    let mut window = Window::new("convex_hull2d").await;
    let mut camera = PanZoomCamera2d::new(Vec2::ZERO, 4.0);
    let mut scene = SceneNode2d::empty();

    let count = 9;
    let mut pts = vec![Vec2::ZERO; count];

    let render_pos = Vec2::ZERO;

    let start_time = web_time::Instant::now();

    while window.render_2d(&mut scene, &mut camera).await {
        let elapsed_time = start_time.elapsed().as_secs_f32();
        let elapsed_time_slow = elapsed_time * 0.2;

        for (i, pt) in pts.iter_mut().enumerate() {
            *pt = lissajous_2d_with_params(
                (i * i) as f32 + elapsed_time_slow,
                2.0 + i as f32 / 3.0,
                (i as f32 / count as f32) + elapsed_time_slow.cos() * 0.1,
                (elapsed_time_slow as f32 + i as f32).cos() * 0.1 + FRAC_PI_2,
                FRAC_PI_4,
            ) * 5f32;
            draw_point(&mut window, *pt, RENDER_SCALE, render_pos, RED);
        }

        /*
         *
         * Compute the convex hull.
         *
         */
        let convex_hull = transformation::convex_hull(&pts);
        draw_polygon(&mut window, &convex_hull, RENDER_SCALE, render_pos, WHITE);
    }
}
