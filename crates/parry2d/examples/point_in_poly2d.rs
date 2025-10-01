mod common_macroquad2d;

use common_macroquad2d::{draw_point, draw_polygon};
use macroquad::prelude::*;
use nalgebra::{Point2, UnitComplex, Vector2};
use parry2d::utils::point_in_poly2d;

const RENDER_SCALE: f32 = 30.0;

#[macroquad::main("points_in_poly2d")]
async fn main() {
    let mut spikes = spikes_polygon();
    let mut squares = squares_polygon();
    let test_points = grid_points();

    let animation_rotation = UnitComplex::new(0.02);
    let polygon_render_pos = Point2::new(screen_width() / 2.0, screen_height() / 2.0);

    for i in 0.. {
        let polygon = if (i / 350) % 2 == 0 {
            &mut spikes
        } else {
            &mut squares
        };

        clear_background(BLACK);

        polygon
            .iter_mut()
            .for_each(|pt| *pt = animation_rotation * *pt);

        draw_polygon(&polygon, RENDER_SCALE, polygon_render_pos, BLUE);

        /*
         * Compute polygon intersections.
         */
        for point in &test_points {
            if point_in_poly2d(point, &polygon) {
                draw_point(*point, RENDER_SCALE, polygon_render_pos, RED);
            } else {
                draw_point(*point, RENDER_SCALE, polygon_render_pos, GREEN);
            }
        }

        next_frame().await
    }
}

fn spikes_polygon() -> Vec<Point2<f32>> {
    let teeths = 3;
    let width = 15.0;
    let height = 7.5;
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
        polygon.push(Point2::new(x + tooth_width / 2.0, height * 1.5) - center);
    }

    polygon
}

fn squares_polygon() -> Vec<Point2<f32>> {
    let scale = 3.0;
    [
        Point2::new(-1.0, -1.0) * scale,
        Point2::new(0.0, -1.0) * scale,
        Point2::new(0.0, 1.0) * scale,
        Point2::new(-2.0, 1.0) * scale,
        Point2::new(-2.0, -2.0) * scale,
        Point2::new(1.0, -2.0) * scale,
        Point2::new(1.0, 2.0) * scale,
        Point2::new(-1.0, 2.0) * scale,
    ]
    .to_vec()
}

fn grid_points() -> Vec<Point2<f32>> {
    let count = 40;
    let spacing = 0.6;
    let mut pts = vec![];
    for i in 0..count {
        for j in 0..count {
            pts.push(Point2::new(
                (i as f32 - count as f32 / 2.0) * spacing,
                (j as f32 - count as f32 / 2.0) * spacing,
            ));
        }
    }
    pts
}
