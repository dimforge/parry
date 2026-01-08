mod utils2d;

use kiss3d::prelude::*;
use parry2d::math::Rotation;
use parry2d::shape::Ball;
use parry2d::transformation::polygons_intersection_points;
use utils2d::draw_polygon;

const RENDER_SCALE: f32 = 30.0;

#[kiss3d::main]
async fn main() {
    let mut window = Window::new("polygons_intersection2d").await;
    let mut camera = PanZoomCamera2d::new(Vec2::ZERO, 2.0);
    let mut scene = SceneNode2d::empty();
    let font = Font::default();

    let spikes = spikes_polygon();
    let mut animated_spikes = spikes.clone();

    let star = star_polygon();

    let animation_scale = 2.0;
    let animation_rotation = Rotation::new(0.008);

    let spikes_render_pos = Vec2::new(-150.0, 0.0);
    let star_render_pos = Vec2::new(150.0, 0.0);

    let start_time = web_time::Instant::now();

    while window.render_2d(&mut scene, &mut camera).await {
        let i = (start_time.elapsed().as_secs_f32() * 60.0) as i32;

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
        draw_polygon(&mut window, &spikes, RENDER_SCALE, spikes_render_pos, BLUE);
        draw_polygon(
            &mut window,
            &animated_spikes,
            RENDER_SCALE,
            spikes_render_pos,
            GREEN,
        );

        draw_polygon(&mut window, &star, RENDER_SCALE, star_render_pos, BLUE);
        draw_polygon(
            &mut window,
            &animated_star,
            RENDER_SCALE,
            star_render_pos,
            GREEN,
        );

        if let Ok(intersections) = spikes_intersections {
            window.draw_text(
                &format!("# spikes intersections: {}", intersections.len()),
                Vec2::new(0.0, 15.0),
                20.0,
                &font,
                WHITE,
            );
            for intersection in intersections {
                draw_polygon(
                    &mut window,
                    &intersection,
                    RENDER_SCALE,
                    spikes_render_pos,
                    RED,
                );
            }
        }

        if let Ok(intersections) = star_intersections {
            window.draw_text(
                &format!("# star intersections: {}", intersections.len()),
                Vec2::new(0.0, 35.0),
                20.0,
                &font,
                WHITE,
            );
            for intersection in intersections {
                draw_polygon(
                    &mut window,
                    &intersection,
                    RENDER_SCALE,
                    star_render_pos,
                    RED,
                );
            }
        }
    }
}

fn star_polygon() -> Vec<Vec2> {
    let mut star = Ball::new(1.5).to_polyline(10);
    star.iter_mut().step_by(2).for_each(|pt| *pt = *pt * 0.6);
    star
}

fn spikes_polygon() -> Vec<Vec2> {
    let teeths = 5;
    let width = 10.0;
    let height = 5.0;
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
        polygon.push(Vec2::new(x + tooth_width / 2.0, height * 0.8) - center);
    }

    polygon
}
