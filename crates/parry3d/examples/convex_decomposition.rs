mod common_macroquad;

extern crate nalgebra as na;

use common_macroquad::mquad_mesh_from_points;
use macroquad::prelude::*;
use na::Vector3;
use parry3d::shape::{Cuboid, SharedShape};

#[macroquad::main("parry2d::utils::point_in_poly2d")]
async fn main() {
    let camera_pos = Vec3::new(8f32, 8f32, 12f32);

    loop {
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
        let cube1 = Cuboid::new(Vector3::repeat(0.5)).to_trimesh();

        dbg!("before convex");
        let convex_mesh = SharedShape::convex_decomposition(&cube1.0, &cube1.1);
        dbg!("after convex, not reached");
        let trimesh_convex = convex_mesh.0.as_trimesh().unwrap();

        /*
         * Display the shapes.
         */
        let trimesh_raw = (
            trimesh_convex.vertices().to_vec(),
            trimesh_convex.indices().to_vec(),
        );
        let mesh = mquad_mesh_from_points(
            &trimesh_raw,
            Vec3::new(1f32, 3f32, 3f32),
            Color::from_rgba(200, 200, 200, 150),
        );
        draw_mesh(&mesh);
        next_frame().await
    }
}
