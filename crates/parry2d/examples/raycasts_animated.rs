mod utils2d;

use kiss3d::prelude::*;
use parry2d::math::{Pose, Rotation};
use parry2d::query::{Ray, RayCast};
use parry2d::shape::Cuboid;
use utils2d::draw_point;

const RENDER_SCALE: f32 = 30.0;

#[kiss3d::main]
async fn main() {
    let mut window = Window::new("raycasts_animated").await;
    let mut camera = PanZoomCamera2d::new(Vec2::ZERO, 4.0);
    let mut scene = SceneNode2d::empty();

    let animation_scale = 1.4;
    let animation_rotation = 0.04;

    let screen_shift = Vec2::ZERO;
    let start_time = web_time::Instant::now();

    while window.render_2d(&mut scene, &mut camera).await {
        let i = (start_time.elapsed().as_secs_f32() * 60.0) as i32 + 1;

        /*
         *
         * Compute the scaled cuboid.
         *
         */
        let cube =
            Cuboid::new(Vec2::new(2.0, 2.0) * ((i as f32 / 50.0).sin().abs() * animation_scale));
        let cube_pose = Pose::new(Vec2::ZERO, 0.008 * i as f32);
        /*
         *
         * Prepare a Raycast and compute its result against the shape.
         *
         */
        let ray = Ray::new(
            Vec2::new(2.0, 2.0),
            Rotation::new(animation_rotation * i as f32).transform_vector(-Vec2::X),
        );
        let toi = cube.cast_ray(&cube_pose, &ray, f32::MAX, true);

        /*
         *
         * Render the raycast's result.
         *
         */
        if let Some(toi) = toi {
            if toi == 0f32 {
                draw_point(&mut window, ray.origin, RENDER_SCALE, screen_shift, YELLOW);
            } else {
                drawline_from_to(
                    &mut window,
                    ray.origin,
                    ray.origin + ray.dir * toi,
                    RENDER_SCALE,
                    screen_shift,
                    GREEN,
                );
            }
        } else {
            drawline_from_to(
                &mut window,
                ray.origin,
                ray.origin + ray.dir * 1000f32,
                RENDER_SCALE,
                screen_shift,
                RED,
            );
        }

        /*
         *
         * Render the cuboid.
         *
         */
        draw_cuboid_polygon(
            &mut window,
            &cube.to_polyline(),
            &cube_pose,
            RENDER_SCALE,
            screen_shift,
            GREEN,
        );
    }
}

fn draw_cuboid_polygon(
    window: &mut Window,
    polygon: &[Vec2],
    pose: &Pose,
    scale: f32,
    shift: Vec2,
    color: Color,
) {
    for i in 0..polygon.len() {
        let a = pose * (polygon[i] * scale) + shift;
        let b = pose * (polygon[(i + 1) % polygon.len()] * scale) + shift;
        window.draw_line_2d(a, b, color, 2.0);
    }
}

fn drawline_from_to(
    window: &mut Window,
    from: Vec2,
    to: Vec2,
    scale: f32,
    shift: Vec2,
    color: Color,
) {
    window.draw_line_2d(from * scale + shift, to * scale + shift, color, 2.0);
}
