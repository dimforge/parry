mod utils3d;

use kiss3d::prelude::*;
use parry3d::bounding_volume::BoundingVolume;
use parry3d::math::Pose;
use parry3d::shape::Ball;
use utils3d::{draw_aabb3, lissajous_3d};

#[kiss3d::main]
async fn main() {
    let mut window = Window::new("aabb3d").await;
    let mut camera = OrbitCamera3d::new(Vec3::new(8.0, 8.0, 12.0), Vec3::new(0.5, 0.0, 0.5));
    let mut scene = SceneNode3d::empty();
    scene
        .add_light(Light::point(100.0))
        .set_position(Vec3::new(0.0, 10.0, 10.0));

    let start_time = web_time::Instant::now();

    while window.render_3d(&mut scene, &mut camera).await {
        let elapsed_time = start_time.elapsed().as_secs_f32() * 0.7;

        /*
         * Initialize the shapes.
         */
        let ball1 = Ball::new(0.5);
        let ball2 = Ball::new(1.0);

        let ball1_translation = lissajous_3d(elapsed_time) * 4f32;
        let ball1_pos = Pose::from_translation(ball1_translation);
        let ball2_pos = Pose::identity();

        /*
         * Compute their axis-aligned bounding boxes.
         */
        let aabb_ball1 = ball1.aabb(&ball1_pos);
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

        window.draw_point(ball1_pos.translation, color, ball1.radius * 50.0);
        window.draw_point(ball2_pos.translation, color, ball2.radius * 50.0);

        draw_aabb3(&mut window, aabb_ball1.mins, aabb_ball1.maxs, color);
        draw_aabb3(&mut window, aabb_ball2.mins, aabb_ball2.maxs, color);
        draw_aabb3(&mut window, bounding_aabb.mins, bounding_aabb.maxs, YELLOW);

        let color_included: Color = if loose_aabb_ball2.contains(&aabb_ball1) {
            BLUE
        } else {
            MAGENTA
        };
        draw_aabb3(
            &mut window,
            loose_aabb_ball2.mins,
            loose_aabb_ball2.maxs,
            color_included,
        );
    }
}
