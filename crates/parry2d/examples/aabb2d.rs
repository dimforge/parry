mod utils2d;

use kiss3d::prelude::*;
use parry2d::bounding_volume::BoundingVolume;
use parry2d::math::Pose;
use parry2d::shape::Ball;
use utils2d::{draw_aabb2, draw_circle, lissajous_2d};

const RENDER_SCALE: f32 = 30.0;

#[kiss3d::main]
async fn main() {
    let mut window = Window::new("aabb2d").await;
    let mut camera = PanZoomCamera2d::new(Vec2::ZERO, 4.0);
    let mut scene = SceneNode2d::empty();

    let start_time = web_time::Instant::now();

    while window.render_2d(&mut scene, &mut camera).await {
        let elapsed_time = start_time.elapsed().as_secs_f32() * 0.7;

        /*
         * Initialize the shapes.
         */
        let ball1 = Ball::new(0.5);
        let ball2 = Ball::new(1.0);

        let ball1_pos = lissajous_2d(elapsed_time) * 5f32;
        let ball1_pose = Pose::from_translation(ball1_pos);
        let ball2_pose = Pose::identity();

        /*
         * Compute their axis-aligned bounding boxes.
         */
        let aabb_ball1 = ball1.aabb(&ball1_pose);
        let aabb_ball2 = ball2.aabb(&ball2_pose);

        // Merge the two boxes.
        let bounding_aabb = aabb_ball1.merged(&aabb_ball2);

        // Enlarge the ball2 aabb.
        let loose_aabb_ball2 = aabb_ball2.loosened(2f32);

        // Intersection test
        let color = if aabb_ball1.intersects(&aabb_ball2) {
            RED
        } else {
            GREEN
        };

        assert!(bounding_aabb.contains(&aabb_ball1));
        assert!(bounding_aabb.contains(&aabb_ball2));
        assert!(loose_aabb_ball2.contains(&aabb_ball2));

        draw_circle(
            &mut window,
            ball1_pos * RENDER_SCALE,
            ball1.radius * RENDER_SCALE,
            color,
        );
        draw_circle(
            &mut window,
            ball2_pose.translation * RENDER_SCALE,
            ball2.radius * RENDER_SCALE,
            color,
        );

        draw_aabb2(
            &mut window,
            aabb_ball1.mins * RENDER_SCALE,
            aabb_ball1.maxs * RENDER_SCALE,
            color,
        );
        draw_aabb2(
            &mut window,
            aabb_ball2.mins * RENDER_SCALE,
            aabb_ball2.maxs * RENDER_SCALE,
            color,
        );
        draw_aabb2(
            &mut window,
            bounding_aabb.mins * RENDER_SCALE,
            bounding_aabb.maxs * RENDER_SCALE,
            YELLOW,
        );

        // Inclusion test
        let color_included: Color = if loose_aabb_ball2.contains(&aabb_ball1) {
            BLUE
        } else {
            MAGENTA
        };
        draw_aabb2(
            &mut window,
            loose_aabb_ball2.mins * RENDER_SCALE,
            loose_aabb_ball2.maxs * RENDER_SCALE,
            color_included,
        );
    }
}
