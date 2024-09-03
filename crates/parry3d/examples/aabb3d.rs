mod common_macroquad;

extern crate nalgebra as na;

use common_macroquad::{draw_polyline, lissajous_3d, mquad_from_na, na_from_mquad};
use macroquad::prelude::*;
use na::Isometry3;
use parry3d::bounding_volume::{Aabb, BoundingVolume};
use parry3d::shape::Ball;
use std::f32::consts::{FRAC_PI_2, FRAC_PI_4};

const RENDER_SCALE: f32 = 30.0;

#[macroquad::main("parry2d::utils::point_in_poly2d")]
async fn main() {
    let render_pos = Vec2::new(300.0, 300.0);
    let camera_pos = Vec3::new(8f32, 8f32, 12f32);

    loop {
        let elapsed_time = get_time() as f32 * 0.7;
        clear_background(BLACK);
        // Initialize 3D camera.
        set_camera(&Camera3D {
            position: camera_pos,
            up: Vec3::new(0f32, 1f32, 0f32),
            target: Vec3::new(0.5f32, 0f32, 0.5f32),
            ..Default::default()
        });

        /*
         * Initialize the shapes.
         */
        let ball1 = Ball::new(0.5);
        let ball2 = Ball::new(1.0);

        let ball1_pos = na_from_mquad(lissajous_3d(elapsed_time)) * 4f32;
        let ball2_pos = Isometry3::identity();

        /*
         * Compute their axis-aligned bounding boxes.
         */
        let aabb_ball1 = ball1.aabb(&ball1_pos.into());
        let aabb_ball2 = ball2.aabb(&ball2_pos);

        // Merge the two boxes.
        let bounding_aabb = aabb_ball1.merged(&aabb_ball2);

        // Enlarge the ball2 aabb.
        let loose_aabb_ball2 = aabb_ball2.loosened(2.25f32);

        // Intersection and inclusion tests.
        let color = if aabb_ball1.intersects(&aabb_ball2) {
            RED
        } else {
            GREEN
        };

        assert!(bounding_aabb.contains(&aabb_ball1));
        assert!(bounding_aabb.contains(&aabb_ball2));
        assert!(loose_aabb_ball2.contains(&aabb_ball2));

        let ball1_translation = mquad_from_na(ball1_pos.coords.into());
        draw_sphere(ball1_translation, ball1.radius, None, color);
        let ball2_translation = mquad_from_na(ball2_pos.translation.vector.into());
        draw_sphere(ball2_translation, ball2.radius, None, color);

        draw_aabb(aabb_ball1, color);
        draw_aabb(aabb_ball2, color);
        draw_aabb(bounding_aabb, YELLOW);

        let color_included: Color = if loose_aabb_ball2.contains(&aabb_ball1) {
            BLUE
        } else {
            MAGENTA
        };
        draw_aabb(loose_aabb_ball2, color_included);
        next_frame().await
    }
}

fn draw_aabb(aabb: Aabb, color: Color) {
    let mins = mquad_from_na(aabb.mins);
    let maxs = mquad_from_na(aabb.maxs);

    let lines = vec![
        (
            Vec3::new(mins.x, mins.y, mins.z), // A -> B
            Vec3::new(maxs.x, mins.y, mins.z),
        ),
        (
            Vec3::new(maxs.x, mins.y, mins.z), // B -> F
            Vec3::new(maxs.x, mins.y, maxs.z),
        ),
        (
            Vec3::new(maxs.x, mins.y, maxs.z), // F -> E
            Vec3::new(mins.x, mins.y, maxs.z),
        ),
        (
            Vec3::new(mins.x, mins.y, maxs.z), // E -> A
            Vec3::new(mins.x, mins.y, mins.z),
        ),
        (
            Vec3::new(mins.x, mins.y, mins.z), // A -> D
            Vec3::new(mins.x, maxs.y, mins.z),
        ),
        (
            Vec3::new(mins.x, maxs.y, mins.z), // D -> H
            Vec3::new(mins.x, maxs.y, maxs.z),
        ),
        (
            Vec3::new(mins.x, maxs.y, maxs.z), // H -> E
            Vec3::new(mins.x, mins.y, maxs.z),
        ),
        (
            Vec3::new(mins.x, maxs.y, maxs.z), // H -> G
            Vec3::new(maxs.x, maxs.y, maxs.z),
        ),
        (
            Vec3::new(maxs.x, maxs.y, maxs.z), // G -> F
            Vec3::new(maxs.x, mins.y, maxs.z),
        ),
        (
            Vec3::new(maxs.x, maxs.y, maxs.z), // G -> C
            Vec3::new(maxs.x, maxs.y, mins.z),
        ),
        (
            Vec3::new(maxs.x, maxs.y, mins.z), // C -> B
            Vec3::new(maxs.x, mins.y, mins.z),
        ),
        (
            Vec3::new(mins.x, maxs.y, mins.z), // D -> C
            Vec3::new(maxs.x, maxs.y, mins.z),
        ),
    ];

    draw_polyline(lines, color);
}
