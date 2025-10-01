mod common_macroquad2d;

extern crate nalgebra as na;

use common_macroquad2d::{draw_polyline, lissajous_2d, mquad_from_na, na_from_mquad};
use macroquad::prelude::*;
use na::{Isometry2, Vector2};
use parry2d::bounding_volume::{Aabb, BoundingVolume};
use parry2d::shape::Cuboid;

const RENDER_SCALE: f32 = 30.0;

#[macroquad::main("bounding_sphere2d")]
async fn main() {
    let render_pos = Vec2::new(300.0, 300.0);

    loop {
        let elapsed_time = get_time() as f32 * 0.7;
        clear_background(BLACK);

        /*
         * Initialize the shapes.
         */
        let cube1: Cuboid = Cuboid::new(Vector2::repeat(0.5));
        let cube2 = Cuboid::new(Vector2::new(1., 0.5));

        let cube1_pos = na_from_mquad(lissajous_2d(elapsed_time)) * 5f32;
        let cube1_pos = Isometry2::from(cube1_pos);
        let cube2_pos = Isometry2::identity();

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

        assert!(loose_bounding_sphere_cube2.contains(&bounding_sphere_cube1));
        assert!(loose_bounding_sphere_cube2.contains(&bounding_sphere_cube2));

        let cube1_translation =
            mquad_from_na(cube1_pos.translation.vector.into()) * RENDER_SCALE + render_pos;
        draw_cuboid(cube1, cube1_translation, color);

        let cube2_translation =
            mquad_from_na(cube2_pos.translation.vector.into()) * RENDER_SCALE + render_pos;
        draw_cuboid(cube2, cube2_translation, color);
        draw_circle_lines(
            bounding_sphere_cube1.center.x * RENDER_SCALE + render_pos.x,
            bounding_sphere_cube1.center.y * RENDER_SCALE + render_pos.y,
            bounding_sphere_cube1.radius * RENDER_SCALE,
            2f32,
            color,
        );
        draw_circle_lines(
            bounding_sphere_cube2.center.x * RENDER_SCALE + render_pos.x,
            bounding_sphere_cube2.center.y * RENDER_SCALE + render_pos.y,
            bounding_sphere_cube2.radius * RENDER_SCALE,
            2f32,
            color,
        );
        draw_circle_lines(
            bounding_bounding_sphere.center.x * RENDER_SCALE + render_pos.x,
            bounding_bounding_sphere.center.y * RENDER_SCALE + render_pos.y,
            bounding_bounding_sphere.radius * RENDER_SCALE,
            2f32,
            YELLOW,
        );

        // Inclusion test
        let color_included: Color = if loose_bounding_sphere_cube2.contains(&bounding_sphere_cube1)
        {
            BLUE
        } else {
            MAGENTA
        };
        draw_circle_lines(
            loose_bounding_sphere_cube2.center.x * RENDER_SCALE + render_pos.x,
            loose_bounding_sphere_cube2.center.y * RENDER_SCALE + render_pos.y,
            loose_bounding_sphere_cube2.radius * RENDER_SCALE,
            2f32,
            color_included,
        );
        next_frame().await
    }
}

fn draw_cuboid(cuboid: Cuboid, pos: Vec2, color: Color) {
    let aabb = cuboid.local_aabb();
    draw_aabb(aabb, pos, color)
}

fn draw_aabb(aabb: Aabb, offset: Vec2, color: Color) {
    let mins = mquad_from_na(aabb.mins) * RENDER_SCALE + offset;
    let maxs = mquad_from_na(aabb.maxs) * RENDER_SCALE + offset;

    let line = vec![
        Vec2::new(mins.x, mins.y),
        Vec2::new(mins.x, maxs.y),
        Vec2::new(maxs.x, maxs.y),
        Vec2::new(maxs.x, mins.y),
        Vec2::new(mins.x, mins.y),
    ];
    let drawable_line = line
        .iter()
        .zip(line.iter().cycle().skip(1).take(line.len()))
        .map(|item| (item.0.clone(), item.1.clone()))
        .collect();
    draw_polyline(drawable_line, color);
}
