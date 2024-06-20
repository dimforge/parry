use macroquad::prelude::*;
use nalgebra::{Point2, UnitComplex, Vector2};
use parry2d::shape::Ball;
use parry2d::transformation::polygons_intersection_points;
use parry2d::utils::point_in_poly2d;

const RENDER_SCALE: f32 = 30.0;

#[macroquad::main("parry2d::utils::point_in_poly2d")]
async fn main() {
    let mut spikes = spikes_polygon();
    let test_points = grid_points();

    let animation_rotation = UnitComplex::new(0.02);
    let spikes_render_pos = Point2::new(screen_width() / 2.0, screen_height() / 2.0);

    for i in 0.. {
        clear_background(BLACK);

        /*
         * Compute polygon intersections.
         */
        spikes
            .iter_mut()
            .for_each(|pt| *pt = animation_rotation * *pt);
        draw_polygon(&spikes, RENDER_SCALE, spikes_render_pos, BLUE);

        for point in &test_points {
            if point_in_poly2d(point, &spikes) {
                draw_point(*point, RENDER_SCALE, spikes_render_pos, RED);
            } else {
                draw_point(*point, RENDER_SCALE, spikes_render_pos, GREEN);
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
