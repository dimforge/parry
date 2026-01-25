mod utils2d;

use kiss3d::prelude::*;
use parry2d::bounding_volume::BoundingVolume;
use parry2d::math::Pose;
use parry2d::shape::Cuboid;
use utils2d::{draw_aabb2, draw_circle, lissajous_2d};

const RENDER_SCALE: f32 = 30.0;

#[kiss3d::main]
async fn main() {
    let mut window = Window::new("bounding_sphere2d").await;
    let mut camera = PanZoomCamera2d::new(Vec2::ZERO, 4.0);
    let mut scene = SceneNode2d::empty();

    let start_time = web_time::Instant::now();

    while window.render_2d(&mut scene, &mut camera).await {
        let elapsed_time = start_time.elapsed().as_secs_f32() * 0.7;

        /*
         * Initialize the shapes.
         */
        let cube1: Cuboid = Cuboid::new(Vec2::splat(0.5));
        let cube2 = Cuboid::new(Vec2::new(1., 0.5));

        let cube1_pt = lissajous_2d(elapsed_time) * 5f32;
        let cube1_pos = Pose::from_translation(cube1_pt);
        let cube2_pos = Pose::identity();

        /*
         * Compute their bounding spheres.
         */
        let bounding_sphere_cube1 = cube1.bounding_sphere(&cube1_pos);
        let bounding_sphere_cube2 = cube2.bounding_sphere(&cube2_pos);

        // Merge the two spheres.
        let bounding_bounding_sphere = bounding_sphere_cube1.merged(&bounding_sphere_cube2);

        // Enlarge the cube2 bounding sphere.
        let loose_bounding_sphere_cube2 = bounding_sphere_cube2.loosened(3.0);

        // Intersection test
        let color = if bounding_sphere_cube1.intersects(&bounding_sphere_cube2) {
            RED
        } else {
            GREEN
        };

        // Due to float imprecisions, it's dangerous to assume that both shapes will be
        // contained in the merged.
        // You can leverage `BoundingVolume::loosened` with an epsilon for expected results.
        //
        // These might fail:
        // assert!(bounding_bounding_sphere.contains(&bounding_sphere_cube1));
        // assert!(bounding_bounding_sphere.contains(&bounding_sphere_cube2));

        // assert!(loose_bounding_sphere_cube2.contains(&bounding_sphere_cube1));
        assert!(loose_bounding_sphere_cube2.contains(&bounding_sphere_cube2));

        // Draw cuboids using their local AABBs
        let aabb1 = cube1.local_aabb();
        let aabb2 = cube2.local_aabb();
        draw_aabb2(
            &mut window,
            aabb1.mins * RENDER_SCALE + cube1_pos.translation * RENDER_SCALE,
            aabb1.maxs * RENDER_SCALE + cube1_pos.translation * RENDER_SCALE,
            color,
        );
        draw_aabb2(
            &mut window,
            aabb2.mins * RENDER_SCALE + cube2_pos.translation * RENDER_SCALE,
            aabb2.maxs * RENDER_SCALE + cube2_pos.translation * RENDER_SCALE,
            color,
        );

        draw_circle(
            &mut window,
            bounding_sphere_cube1.center * RENDER_SCALE,
            bounding_sphere_cube1.radius * RENDER_SCALE,
            color,
        );
        draw_circle(
            &mut window,
            bounding_sphere_cube2.center * RENDER_SCALE,
            bounding_sphere_cube2.radius * RENDER_SCALE,
            color,
        );
        draw_circle(
            &mut window,
            bounding_bounding_sphere.center * RENDER_SCALE,
            bounding_bounding_sphere.radius * RENDER_SCALE,
            YELLOW,
        );

        // Inclusion test
        let color_included: Color = if loose_bounding_sphere_cube2.contains(&bounding_sphere_cube1)
        {
            BLUE
        } else {
            MAGENTA
        };
        draw_circle(
            &mut window,
            loose_bounding_sphere_cube2.center * RENDER_SCALE,
            loose_bounding_sphere_cube2.radius * RENDER_SCALE,
            color_included,
        );
    }
}
