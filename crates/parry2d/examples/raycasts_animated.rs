use macroquad::prelude::*;
use nalgebra::{Isometry2, Point2, UnitComplex, Vector2};
use parry2d::math::Isometry;
use parry2d::query::{Ray, RayCast};
use parry2d::shape::Cuboid;

const RENDER_SCALE: f32 = 30.0;

#[macroquad::main("parry2d::query::RayCast")]
async fn main() {
    let animation_scale = 1.4;
    let animation_rotation = 0.04;

    for i in 1.. {
        clear_background(BLACK);

        let screen_shift = Point2::new(screen_width() / 2.0, screen_height() / 2.0);
        /*
         *
         * Compute the scaled cuboid.
         *
         */
        let cube =
            Cuboid::new(Vector2::new(2.0, 2.0) * ((i as f32 / 50.0).sin().abs() * animation_scale));
        let cube_pose = Isometry2::rotation(0.008 * i as f32);
        /*
         *
         * Prepare a Raycast and compute its result against the shape.
         *
         */
        let ray = Ray::new(
            Point2::new(2.0, 2.0),
            UnitComplex::new(animation_rotation * i as f32) * -Vector2::x(),
        );
        let toi = cube.cast_ray(&cube_pose, &ray, std::f32::MAX, true);

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

fn draw_polygon(
    polygon: &[Point2<f32>],
    pose: &Isometry<f32>,
    scale: f32,
    shift: Point2<f32>,
    color: Color,
) {
    for i in 0..polygon.len() {
        let a = pose * (scale * polygon[i]);
        let b = pose * (scale * polygon[(i + 1) % polygon.len()]);
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

fn draw_point(point: Point2<f32>, scale: f32, shift: Point2<f32>, color: Color) {
    let edge_len = 0.15;
    draw_line(
        (point.x - edge_len) * scale + shift.x,
        point.y * scale + shift.y,
        (point.x + edge_len) * scale + shift.x,
        point.y * scale + shift.y,
        2.0,
        color,
    );
    draw_line(
        point.x * scale + shift.x,
        (point.y - edge_len) * scale + shift.y,
        point.x * scale + shift.x,
        (point.y + edge_len) * scale + shift.y,
        2.0,
        color,
    );
}

fn drawline_from_to(
    from: Point2<f32>,
    to: Point2<f32>,
    scale: f32,
    shift: Point2<f32>,
    color: Color,
) {
    let from = from * scale + shift.coords;
    let to = to * scale + shift.coords;
    draw_line(from.x, from.y, to.x, to.y, 2.0, color);
}
