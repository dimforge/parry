mod common_macroquad3d;

extern crate nalgebra as na;

use std::ops::Rem;

use common_macroquad3d::{lissajous_3d, mquad_from_na, na_from_mquad};
use macroquad::prelude::*;
use na::{Isometry3, Vector3};
use parry3d::bounding_volume::BoundingVolume;
use parry3d::shape::Cuboid;

#[macroquad::main("bounding_sphere3d")]
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
        let cube1 = Cuboid::new(Vector3::repeat(0.5));
        let cube2 = Cuboid::new(Vector3::new(0.5, 1.0, 0.5));

        let cube1_pos = na_from_mquad(lissajous_3d(elapsed_time)) * 4f32;
        let cube1_pos = Isometry3::from(cube1_pos);
        let cube2_pos = Isometry3::identity(); // Identity matrix.

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
        let mut color = if bounding_sphere_cube1.intersects(&bounding_sphere_cube2) {
            RED
        } else {
            GREEN
        };
        color.a = 1f32 * (elapsed_time.rem(1f32) - 0.5).abs() * 2f32;

        // Due to float imprecisions, it's dangerous to assume that both shapes will be
        // contained in the merged.
        // You can leverage `BoundingVolume::loosened` with an epsilon for expected results.
        //
        // These might fail:
        //assert!(bounding_bounding_sphere.contains(&bounding_sphere_cube1));
        //assert!(bounding_bounding_sphere.contains(&bounding_sphere_cube2));
        assert!(loose_bounding_sphere_cube2.contains(&bounding_sphere_cube2));

        let cube1_translation = mquad_from_na(cube1_pos.translation.vector.into());
        draw_cube_wires(
            cube1_translation,
            mquad_from_na(cube1.half_extents.into()) * 2f32,
            WHITE,
        );
        let cube2_translation = mquad_from_na(cube2_pos.translation.vector.into());
        draw_cube_wires(
            cube2_translation,
            mquad_from_na(cube2.half_extents.into()) * 2f32,
            WHITE,
        );

        draw_sphere_wires(
            mquad_from_na(bounding_sphere_cube1.center),
            bounding_sphere_cube1.radius,
            None,
            color,
        );
        draw_sphere_wires(
            mquad_from_na(bounding_sphere_cube2.center),
            bounding_sphere_cube2.radius,
            None,
            color,
        );
        draw_sphere_wires(
            mquad_from_na(bounding_bounding_sphere.center),
            bounding_bounding_sphere.radius,
            None,
            YELLOW,
        );

        let color_included: Color = if loose_bounding_sphere_cube2.contains(&bounding_sphere_cube1)
        {
            BLUE
        } else {
            MAGENTA
        };
        draw_sphere_wires(
            mquad_from_na(loose_bounding_sphere_cube2.center),
            loose_bounding_sphere_cube2.radius,
            None,
            color_included,
        );
        next_frame().await
    }
}
