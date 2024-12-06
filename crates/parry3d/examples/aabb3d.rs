mod common_macroquad3d;

extern crate nalgebra as na;

use common_macroquad3d::{lissajous_3d, mquad_from_na, na_from_mquad};
use macroquad::prelude::*;
use na::Isometry3;
use parry3d::bounding_volume::{Aabb, BoundingVolume};
use parry3d::shape::Ball;

#[macroquad::main("aabb3d")]
async fn main() {
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
    let size = aabb.maxs - aabb.mins;
    draw_cube_wires(
        mquad_from_na(aabb.maxs - size / 2f32),
        mquad_from_na(size.into()),
        color,
    );
}
