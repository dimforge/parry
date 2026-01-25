mod utils2d;

use kiss3d::prelude::*;
use parry2d::math::Rotation;
use parry2d::utils::point_in_poly2d;
use utils2d::{draw_point, draw_polygon};

const RENDER_SCALE: f32 = 30.0;

#[kiss3d::main]
async fn main() {
    let mut window = Window::new("points_in_poly2d").await;
    let mut camera = PanZoomCamera2d::new(Vec2::ZERO, 2.0);
    let mut scene = SceneNode2d::empty();

    let mut spikes = spikes_polygon();
    let mut squares = squares_polygon();
    let test_points = grid_points();

    let animation_rotation = Rotation::new(0.02);
    let polygon_render_pos = Vec2::ZERO;

    let start_time = web_time::Instant::now();

    while window.render_2d(&mut scene, &mut camera).await {
        let i = (start_time.elapsed().as_secs_f32() * 60.0) as i32;

        let polygon = if (i / 350) % 2 == 0 {
            &mut spikes
        } else {
            &mut squares
        };

        polygon
            .iter_mut()
            .for_each(|pt| *pt = animation_rotation.transform_vector(*pt));

        draw_polygon(
            &mut window,
            &polygon,
            RENDER_SCALE,
            polygon_render_pos,
            BLUE,
        );

        /*
         * Compute polygon intersections.
         */
        for point in &test_points {
            if point_in_poly2d(*point, &polygon) {
                draw_point(&mut window, *point, RENDER_SCALE, polygon_render_pos, RED);
            } else {
                draw_point(&mut window, *point, RENDER_SCALE, polygon_render_pos, GREEN);
            }
        }
    }
}

fn spikes_polygon() -> Vec<Vec2> {
    let teeths = 3;
    let width = 15.0;
    let height = 7.5;
    let tooth_width = width / (teeths as f32);
    let center = Vec2::new(width / 2.0, height / 2.0);

    let mut polygon = vec![
        Vec2::new(width, 0.0) - center,
        Vec2::new(width, height) - center,
        Vec2::new(0.0, height) - center,
    ];

    for i in 0..teeths {
        let x = i as f32 * tooth_width;
        polygon.push(Vec2::new(x, 0.0) - center);
        polygon.push(Vec2::new(x + tooth_width / 2.0, height * 1.5) - center);
    }

    polygon
}

fn squares_polygon() -> Vec<Vec2> {
    let scale = 3.0;
    [
        Vec2::new(-1.0, -1.0) * scale,
        Vec2::new(0.0, -1.0) * scale,
        Vec2::new(0.0, 1.0) * scale,
        Vec2::new(-2.0, 1.0) * scale,
        Vec2::new(-2.0, -2.0) * scale,
        Vec2::new(1.0, -2.0) * scale,
        Vec2::new(1.0, 2.0) * scale,
        Vec2::new(-1.0, 2.0) * scale,
    ]
    .to_vec()
}

fn grid_points() -> Vec<Vec2> {
    let count = 40;
    let spacing = 0.6;
    let mut pts = vec![];
    for i in 0..count {
        for j in 0..count {
            pts.push(Vec2::new(
                (i as f32 - count as f32 / 2.0) * spacing,
                (j as f32 - count as f32 / 2.0) * spacing,
            ));
        }
    }
    pts
}
