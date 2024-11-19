mod common_macroquad3d;

extern crate nalgebra as na;

use common_macroquad3d::mquad_mesh_from_points;
use macroquad::prelude::*;
use parry3d::{math::Point, shape::SharedShape};

#[macroquad::main("parry2d::utils::point_in_poly2d")]
async fn main() {
    let camera_pos = Vec3::new(72f32, 72f32, 18f32);

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

        //// This cube1 ends up in infinite loop if attempting to call `convex_decomposition` on it
        //
        // let cube1 = parry3d::shape::Cuboid::new(na::Vector3::repeat(0.5)).to_trimesh();
        // dbg!("before convex");
        // let convex_mesh = SharedShape::convex_decomposition(&cube1.0, &cube1.1);
        // dbg!("after convex, not reached");
        // let trimesh_convex = convex_mesh.0.as_trimesh().unwrap();
        //
        ////

        //
        // Those values are from https://github.com/dimforge/rapier/issues/223#issuecomment-927816118

        let points = [
            Point::new(-0.74, -1.74, 52.3025),
            Point::new(-0.74, -1.74, -52.3025),
            Point::new(0.74, -1.74, -52.3025),
            Point::new(0.74, -1.74, 52.3025),
            Point::new(-0.74, 1.74, 52.3025),
            Point::new(-0.74, 1.74, -52.3025),
            Point::new(0.74, 1.74, -52.3025),
            Point::new(0.74, 1.74, 52.3025),
        ];
        let indices = [
            [4, 5, 0],
            [5, 1, 0],
            [5, 6, 1],
            [6, 2, 1],
            [6, 7, 3],
            [2, 6, 3],
            [7, 4, 0],
            [3, 7, 0],
            [0, 1, 2],
            [3, 0, 2],
            [7, 6, 5],
            [4, 7, 5],
        ];

        let convex_mesh = SharedShape::convex_decomposition(&points, &indices);
        let trimesh_convex_compound = convex_mesh.as_compound().unwrap();

        for s in trimesh_convex_compound.shapes() {
            let trimesh_convex = s.1.as_convex_polyhedron().unwrap().to_trimesh();

            /*
             * Display the shapes.
             */
            let mesh = mquad_mesh_from_points(
                &trimesh_convex,
                Vec3::new(1f32, 3f32, 3f32),
                Color::from_rgba(200, 200, 200, 150),
            );
            draw_mesh(&mesh);
        }
        next_frame().await
    }
}
