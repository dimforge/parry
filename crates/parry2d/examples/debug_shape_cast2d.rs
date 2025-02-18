mod common_macroquad2d;

use common_macroquad2d::draw_point;
use macroquad::prelude::*;
use nalgebra::{Isometry2, Point2};
use parry2d::math::Isometry;
use parry2d::query::{self, Ray, ShapeCastOptions};
use parry2d::shape::{Ball, ConvexPolygon, Shape};

const RENDER_SCALE: f32 = 1.0;
const BALLCAST_WIDTH: f32 = 16.0;

#[macroquad::main("raycasts_animated")]
async fn main() {
    for _i in 1.. {
        clear_background(BLACK);

        let screen_shift = Point2::new(screen_width() / 2.0, screen_height() / 2.0);

        let to_cast_against = ConvexPolygon::from_convex_polyline(
            [
                [-24.0, 0.0].into(),
                [0.0, 24.0].into(),
                [24.0, 0.0].into(),
                [0.0, -24.0].into(),
            ]
            .into(),
        )
        .expect("Failed to create ConvexPolygon from polyline");
        let to_cast_against_pose = Isometry2::rotation(0f32);

        let mouse_pos = mouse_position();
        let mouse_position_world =
            (Point2::<f32>::new(mouse_pos.0, mouse_pos.1) - screen_shift.coords) / RENDER_SCALE;
        let target_pos: Point2<f32> = [-312.0, 152.0].into();

        // Those 2 fail with `min_bound >= _eps_tol`, fixed with a tolerance * 100
        shape_cast_debug(
            screen_shift,
            [99.0, -33.0].into(),
            target_pos,
            to_cast_against.clone(),
        );
        shape_cast_debug(
            screen_shift,
            [98.0, -31.0].into(),
            target_pos,
            to_cast_against.clone(),
        );
        // This fails with `niter == 100` (and `niter == 100_000`), fixed with a tolerance * 10_000
        shape_cast_debug(
            screen_shift,
            [47.0, -32.0].into(),
            target_pos,
            to_cast_against.clone(),
        );

        // For debug purposes, raycast to mouse position.
        // Rendered last to be on top of the other raycasts
        shape_cast_debug(
            screen_shift,
            target_pos,
            mouse_position_world,
            to_cast_against.clone(),
        );

        /*
         *
         * Render the cuboid.
         *
         */
        draw_polygon(
            &to_cast_against.points(),
            &to_cast_against_pose,
            RENDER_SCALE,
            screen_shift,
            GREEN,
        );

        next_frame().await
    }
}

fn shape_cast_debug(
    screen_shift: Point2<f32>,
    source_pos: Point2<f32>,
    target_pos: Point2<f32>,
    to_cast_against: ConvexPolygon,
) {
    /*
     *
     * Prepare a Raycast and compute its result against the shape.
     *
     */
    let ray = Ray::new(source_pos, target_pos - source_pos);

    let pos1: Point2<f32> = source_pos.into();
    let vel1 = target_pos - source_pos;
    let g1 = Ball::new(BALLCAST_WIDTH);
    let pos2 = [0.0, 0.0];
    let vel2 = [0.0, 0.0];
    let g2 = to_cast_against.clone_dyn();

    let toi = query::cast_shapes(
        &pos1.into(),
        &vel1,
        &g1,
        &pos2.into(),
        &vel2.into(),
        &*g2,
        ShapeCastOptions::with_max_time_of_impact(1.0),
    )
    .unwrap();

    /*
     *
     * Render the raycast's result.
     *
     */
    drawcircle_at(pos1, BALLCAST_WIDTH, RENDER_SCALE, screen_shift, ORANGE);

    if let Some(toi) = toi {
        if toi.time_of_impact == 0f32 {
            draw_point(ray.origin, RENDER_SCALE, screen_shift, YELLOW);
        } else {
            drawline_from_to(
                ray.origin,
                (ray.point_at(toi.time_of_impact).coords).into(),
                RENDER_SCALE,
                screen_shift,
                GREEN,
            );
        }
        drawcircle_at(
            (ray.point_at(toi.time_of_impact).coords).into(),
            BALLCAST_WIDTH,
            RENDER_SCALE,
            screen_shift,
            GREEN,
        );
    } else {
        drawline_from_to(
            ray.origin,
            ray.origin + ray.dir,
            RENDER_SCALE,
            screen_shift,
            RED,
        );
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

fn drawcircle_at(center: Point2<f32>, radius: f32, scale: f32, shift: Point2<f32>, color: Color) {
    let center = center * scale + shift.coords;
    draw_circle_lines(center.x, center.y, radius, 1f32, color);
}
