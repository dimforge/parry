mod utils3d;

use kiss3d::prelude::*;
use parry3d::bounding_volume::BoundingVolume;
use parry3d::math::Pose;
use parry3d::shape::Cuboid;
use utils3d::{draw_aabb3, draw_sphere_wires, lissajous_3d};

#[kiss3d::main]
async fn main() {
    let mut window = Window::new("bounding_sphere3d").await;
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
        let cube1 = Cuboid::new(Vec3::splat(0.5));
        let cube2 = Cuboid::new(Vec3::new(0.5, 1.0, 0.5));

        let cube1_translation = lissajous_3d(elapsed_time) * 4f32;
        let cube1_pos = Pose::from_translation(cube1_translation);
        let cube2_pos = Pose::identity(); // Identity matrix.

        /*
         * Compute their bounding spheres.
         */
        let bounding_sphere_cube1 = cube1.bounding_sphere(&cube1_pos);
        let bounding_sphere_cube2 = cube2.bounding_sphere(&cube2_pos);

        // Merge the two spheres.
        let bounding_bounding_sphere = bounding_sphere_cube1.merged(&bounding_sphere_cube2);

        // Enlarge the cube2 bounding sphere.
        let loose_bounding_sphere_cube2 = bounding_sphere_cube2.loosened(3.0);

        // Intersection and inclusion tests.
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
        //assert!(bounding_bounding_sphere.contains(&bounding_sphere_cube1));
        //assert!(bounding_bounding_sphere.contains(&bounding_sphere_cube2));
        assert!(loose_bounding_sphere_cube2.contains(&bounding_sphere_cube2));

        // Draw cube wireframes
        let mins1 = cube1_pos.translation - cube1.half_extents;
        let maxs1 = cube1_pos.translation + cube1.half_extents;
        draw_aabb3(&mut window, mins1, maxs1, WHITE);

        let mins2 = cube2_pos.translation - cube2.half_extents;
        let maxs2 = cube2_pos.translation + cube2.half_extents;
        draw_aabb3(&mut window, mins2, maxs2, WHITE);

        draw_sphere_wires(
            &mut window,
            bounding_sphere_cube1.center,
            bounding_sphere_cube1.radius,
            color,
        );
        draw_sphere_wires(
            &mut window,
            bounding_sphere_cube2.center,
            bounding_sphere_cube2.radius,
            color,
        );
        draw_sphere_wires(
            &mut window,
            bounding_bounding_sphere.center,
            bounding_bounding_sphere.radius,
            YELLOW,
        );

        let color_included: Color = if loose_bounding_sphere_cube2.contains(&bounding_sphere_cube1)
        {
            BLUE
        } else {
            MAGENTA
        };
        draw_sphere_wires(
            &mut window,
            loose_bounding_sphere_cube2.center,
            loose_bounding_sphere_cube2.radius,
            color_included,
        );
    }
}
