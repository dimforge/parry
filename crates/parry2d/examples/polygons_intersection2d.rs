mod common_macroquad2d;

use common_macroquad2d::draw_polygon;
use macroquad::prelude::*;
use parry2d::math::{Rotation, Vector, Vector};
use parry2d::shape::Ball;
use parry2d::transformation::polygons_intersection_points;

const RENDER_SCALE: f32 = 30.0;

#[macroquad::main("polygons_intersection2d")]
async fn main() {
    let spikes = spikes_polygon();
    let mut animated_spikes = spikes.clone();

    let star = star_polygon();

    let animation_scale = 2.0;
    let animation_rotation = Rotation::new(0.008);

    let spikes_render_pos = Vector::new(300.0, 300.0);
    let star_render_pos = Vector::new(600.0, 300.0);

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
            .for_each(|pt| *pt = animation_rotation.transform_vector(*pt));
        let spikes_intersections = polygons_intersection_points(&spikes, &animated_spikes);

        let animated_star: Vec<_> = star
            .iter()
            .map(|pt| {
                animation_rotation.powf(i as f32).transform_vector(*pt)
                    * ((i as f32 / 100.0).sin().abs() * animation_scale)
            })
            .collect();

        let star_intersections = polygons_intersection_points(&star, &animated_star);

        /*
         *
         * Render the polygons and their intersections.
         *
         */
        draw_polygon(&spikes, RENDER_SCALE, spikes_render_pos, BLUE);
        draw_polygon(&animated_spikes, RENDER_SCALE, spikes_render_pos, GREEN);

        draw_polygon(&star, RENDER_SCALE, star_render_pos, BLUE);
        draw_polygon(&animated_star, RENDER_SCALE, star_render_pos, GREEN);

        if let Ok(intersections) = spikes_intersections {
            draw_text(
                &format!("# spikes intersections: {}", intersections.len()),
                0.0,
                15.0,
                20.0,
                WHITE,
            );
            for intersection in intersections {
                draw_polygon(&intersection, RENDER_SCALE, spikes_render_pos, RED);
            }
        }

        if let Ok(intersections) = star_intersections {
            draw_text(
                &format!("# star intersections: {}", intersections.len()),
                0.0,
                30.0,
                20.0,
                WHITE,
            );
            for intersection in intersections {
                draw_polygon(&intersection, RENDER_SCALE, star_render_pos, RED);
            }
        }

        next_frame().await
    }
}

fn star_polygon() -> Vec<Vector> {
    let mut star = Ball::new(1.5).to_polyline(10);
    star.iter_mut().step_by(2).for_each(|pt| *pt = *pt * 0.6);
    star
}

fn spikes_polygon() -> Vec<Vector> {
    let teeths = 5;
    let width = 10.0;
    let height = 5.0;
    let tooth_width = width / (teeths as f32);
    let center = Vector::new(width / 2.0, height / 2.0);

    let mut polygon = vec![
        Vector::new(width, 0.0) - center,
        Vector::new(width, height) - center,
        Vector::new(0.0, height) - center,
    ];

    for i in 0..teeths {
        let x = i as f32 * tooth_width;
        polygon.push(Vector::new(x, 0.0) - center);
        polygon.push(Vector::new(x + tooth_width / 2.0, height * 0.8) - center);
    }

    polygon
}
