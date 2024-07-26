use macroquad::prelude::*;
use nalgebra::{Isometry2, Point2, UnitComplex, Vector2};
use parry2d::query::{Ray, RayCast};
use parry2d::shape::{Cuboid, TriMesh};

const RENDER_SCALE: f32 = 30.0;

#[macroquad::main("parry2d::query::RayCast")]
async fn main() {
    let animation_scale = 2.0;
    let animation_rotation = UnitComplex::new(0.008);
    let spikes_render_pos = Point2::new(300.0, 300.0);

    for i in 1.. {
        clear_background(BLACK);

        /*
         *
         * Compute the rotated/scaled cuboid.
         *
         */
        let mut animated_cube = Cuboid::new(Vector2::new(2.0, 2.0)).to_polyline();
        animated_cube.iter_mut().for_each(|pt| {
            *pt = animation_rotation.powf(i as f32)
                * *pt
                * ((i as f32 / 100.0).sin().abs() * animation_scale);
        });
        let trimesh = TriMesh::from_polygon(animated_cube).unwrap();
        /*
         *
         * Prepare a Raycast and compute its result against the shape.
         *
         */
        let ray = Ray::new(
            Point2::new(2.0, 1.0),
            animation_rotation.powf(i as f32 * 2f32) * Vector2::y(),
        );
        let toi = trimesh.cast_ray(&Isometry2::identity(), &ray, std::f32::MAX, true);

        /*
         *
         * Render the raycast's result.
         *
         */
        if let Some(toi) = toi {
            if toi == 0f32 {
                draw_point(ray.origin, RENDER_SCALE, spikes_render_pos, YELLOW);
            } else {
                drawline_from_to(
                    ray.origin,
                    ray.origin + ray.dir * toi,
                    RENDER_SCALE,
                    spikes_render_pos,
                    GREEN,
                );
            }
        } else {
            drawline_from_to(
                ray.origin,
                ray.origin + ray.dir * 1000f32,
                RENDER_SCALE,
                spikes_render_pos,
                RED,
            );
        }

        /*
         *
         * Render the polygons and their intersections.
         *
         */
        draw_polygon(&trimesh.vertices(), RENDER_SCALE, spikes_render_pos, GREEN);

        next_frame().await
    }
}

fn draw_polygon(polygon: &[Point2<f32>], scale: f32, shift: Point2<f32>, color: Color) {
    for i in 0..polygon.len() {
        let a = polygon[i];
        let b = polygon[(i + 1) % polygon.len()];
        draw_line(
            a.x * scale + shift.x,
            a.y * scale + shift.y,
            b.x * scale + shift.x,
            b.y * scale + shift.y,
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
