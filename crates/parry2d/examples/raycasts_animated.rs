mod common_macroquad2d;

use common_macroquad2d::draw_point;
use macroquad::prelude::*;
use parry2d::math::{Pose, Rotation, Vector, Vector};
use parry2d::query::{Ray, RayCast};
use parry2d::shape::Cuboid;

const RENDER_SCALE: f32 = 30.0;

#[macroquad::main("raycasts_animated")]
async fn main() {
    let animation_scale = 1.4;
    let animation_rotation = 0.04;

    for i in 1.. {
        clear_background(BLACK);

        let screen_shift = Vector::new(screen_width() / 2.0, screen_height() / 2.0);
        /*
         *
         * Compute the scaled cuboid.
         *
         */
        let cube =
            Cuboid::new(Vector::new(2.0, 2.0) * ((i as f32 / 50.0).sin().abs() * animation_scale));
        let cube_pose = Pose::new(Vector::ZERO, 0.008 * i as f32);
        /*
         *
         * Prepare a Raycast and compute its result against the shape.
         *
         */
        let ray = Ray::new(
            Vector::new(2.0, 2.0),
            Rotation::new(animation_rotation * i as f32).transform_vector(-Vector::X),
        );
        let toi = cube.cast_ray(&cube_pose, &ray, f32::MAX, true);

        /*
         *
         * Render the raycast's result.
         *
         */
        if let Some(toi) = toi {
            if toi == 0f32 {
                draw_point(ray.origin, RENDER_SCALE, screen_shift, YELLOW);
            } else {
                drawline_from_to(
                    ray.origin,
                    ray.origin + ray.dir * toi,
                    RENDER_SCALE,
                    screen_shift,
                    GREEN,
                );
            }
        } else {
            drawline_from_to(
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
        draw_polygon(
            &cube.to_polyline(),
            &cube_pose,
            RENDER_SCALE,
            screen_shift,
            GREEN,
        );

        next_frame().await
    }
}

fn draw_polygon(polygon: &[Vector], pose: &Pose, scale: f32, shift: Vector, color: Color) {
    for i in 0..polygon.len() {
        let a = pose * (polygon[i] * scale);
        let b = pose * (polygon[(i + 1) % polygon.len()] * scale);
        draw_line(
            a.x + shift.x,
            a.y + shift.y,
            b.x + shift.x,
            b.y + shift.y,
            2.0,
            color,
        );
    }
}

fn drawline_from_to(from: Vector, to: Vector, scale: f32, shift: Vector, color: Color) {
    let from = from * scale + shift;
    let to = to * scale + shift;
    draw_line(from.x, from.y, to.x, to.y, 2.0, color);
}
