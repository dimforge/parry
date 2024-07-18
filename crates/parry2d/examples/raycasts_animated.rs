use macroquad::prelude::*;
use nalgebra::{Isometry2, Point2, UnitComplex, Vector2};
use parry2d::query::{Ray, RayCast};
use parry2d::shape::{Ball, TriMesh};
use parry2d::transformation::polygons_intersection_points;

const RENDER_SCALE: f32 = 30.0;

#[macroquad::main("parry2d::utils::polygons_intersection_points")]
async fn main() {
    let mut animated_spikes = spikes_polygon();

    let mut animated_star = star_polygon();

    let animation_scale = 2.0;
    let animation_rotation = UnitComplex::new(0.008);

    let spikes_render_pos = Point2::new(300.0, 300.0);
    let star_render_pos = Point2::new(600.0, 300.0);

    for i in 0.. {
        clear_background(BLACK);

        /*
         *
         * Compute the rotated/scaled polygons, and compute their intersection with the original
         * polygon.
         *
         */
        animated_spikes
            .iter_mut()
            .for_each(|pt| *pt = animation_rotation * *pt);

        animated_star.iter_mut().for_each(|pt| {
            *pt = animation_rotation.powf(i as f32)
                * *pt
                * ((i as f32 / 100.0).sin().abs() * animation_scale)
        });

        let trimesh = TriMesh::from_polygon(animated_spikes.clone()).unwrap();
        let ray = Ray::new(
            spikes_render_pos.into(),
            animation_rotation.powf(-i as f32) * Vector2::y(),
        );
        let toi = trimesh.cast_ray(&Isometry2::identity(), &ray, std::f32::MAX, false);

        if let Some(toi) = toi {
            drawline_from_to(ray.origin, ray.origin + ray.dir * toi, GREEN);
        } else {
            drawline_from_to(ray.origin, ray.origin + ray.dir * 1000f32, RED);
        }

        /*
         *
         * Render the polygons and their intersections.
         *
         */
        draw_polygon(&trimesh.vertices(), RENDER_SCALE, spikes_render_pos, GREEN);

        draw_polygon(&animated_star, RENDER_SCALE, star_render_pos, GREEN);

        next_frame().await
    }
}

fn star_polygon() -> Vec<Point2<f32>> {
    let mut star = Ball::new(1.5).to_polyline(10);
    star.iter_mut().step_by(2).for_each(|pt| *pt = *pt * 0.6);
    star
}

fn spikes_polygon() -> Vec<Point2<f32>> {
    let teeths = 5;
    let width = 10.0;
    let height = 5.0;
    let tooth_width = width / (teeths as f32);
    let center = Vector2::new(width / 2.0, height / 2.0);

    let mut polygon = vec![
        Point2::new(width, 0.0) - center,
        Point2::new(width, height) - center,
        Point2::new(0.0, height) - center,
    ];

    for i in 0..teeths {
        let x = i as f32 * tooth_width;
        polygon.push(Point2::new(x, 0.0) - center);
        polygon.push(Point2::new(x + tooth_width / 2.0, height * 0.8) - center);
    }

    polygon
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

fn drawline_from_to(from: Point2<f32>, to: Point2<f32>, color: Color) {
    draw_line(from.x, from.y, to.x, to.y, 2.0, color);
}
