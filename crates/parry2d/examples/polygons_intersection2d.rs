use macroquad::prelude::*;
use nalgebra::{Point2, UnitComplex, Vector2};
use parry2d::transformation::polygons_intersection_points;

#[macroquad::main("parry2d::utils::polygons_intersection_points")]
async fn main() {
    let teeths = 10;
    let width = 10.0;
    let height = 5.0;
    let tooth_width = width / (teeths as f32);
    let center = Vector2::new(width / 2.0, height / 2.0);

    let mut polygon = vec![
        Point2::new(width, 0.0) - center,
        Point2::new(width, height) - center,
        Point2::new(0.0, height) - center,
    ];

    let teeths = 5;
    for i in 0..teeths {
        let x = i as f32 * tooth_width;
        polygon.push(Point2::new(x, 0.0) - center);
        polygon.push(Point2::new(x + tooth_width / 2.0, height * 0.8) - center);
    }

    const RENDER_SCALE: f32 = 30.0;
    let mut rotated_polygon = polygon.clone();
    let rot = UnitComplex::new(0.008);
    let shift = Point2::new(300.0, 300.0);

    for i in 0.. {
        println!("====================");
        clear_background(BLACK);

        rotated_polygon.iter_mut().for_each(|pt| *pt = rot * *pt);

        draw_polygon(&polygon, RENDER_SCALE, shift, BLUE);
        draw_polygon(&rotated_polygon, RENDER_SCALE, shift, GREEN);

        if let Ok(intersections) = polygons_intersection_points(&polygon, &rotated_polygon) {
            println!("Found num intersections: {}", intersections.len());
            for intersection in intersections {
                println!("Num vertices: {}", intersection.len());
                draw_polygon(&intersection, RENDER_SCALE, shift, RED);
            }
        } else {
            eprintln!("Entered inifinite loop.");
        }

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
