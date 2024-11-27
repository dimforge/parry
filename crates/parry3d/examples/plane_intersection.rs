use macroquad::prelude::*;
use nalgebra::{UnitVector3, Vector3};
use parry3d::query::IntersectResult;
use parry3d::shape::{Cuboid, TriMesh};

mod common_macroquad3d;
use common_macroquad3d::*;

#[macroquad::main("plane_intersection")]
async fn main() {
    let trimesh = Cuboid::new(Vector3::repeat(1.0)).to_trimesh();

    let light_pos = Vec3::new(-1f32, 3.5f32, -3f32);
    let camera_pos = Vec3::new(-1.5f32, 2.5f32, -3f32);

    let mesh = mquad_mesh_from_points(&trimesh, light_pos, DARKGRAY);
    let trimesh = TriMesh::new(trimesh.0, trimesh.1).unwrap();

    for _ in 1.. {
        clear_background(BLACK);

        let elapsed_time = get_time();

        // Animated rotation for the intersection plane.
        let bias = -1.2 * (elapsed_time as f32 / 3f32).sin();
        let rotation = Quat::from_axis_angle(Vec3::Z, (elapsed_time as f32 * 40f32).to_radians());
        let up_plane_vector = rotation * Vec3::Y;

        // Get the intersection polyline.
        let intersection_result = trimesh.intersection_with_local_plane(
            &UnitVector3::new_normalize(Vector3::new(
                up_plane_vector.x,
                up_plane_vector.y,
                up_plane_vector.z,
            )),
            bias,
            0.0005,
        );

        // Initialize 3D camera.
        set_camera(&Camera3D {
            position: camera_pos,
            up: Vec3::new(0f32, 1f32, 0f32),
            target: Vec3::new(0.5f32, 0f32, 0.5f32),
            ..Default::default()
        });

        // Draw involved shapes.
        let plane_center = up_plane_vector * bias;
        draw_line_3d(plane_center, plane_center + up_plane_vector, GREEN);
        draw_mesh(&mesh);
        draw_grid_ex(10, 0.333, BLUE, RED, plane_center, rotation);

        /*
         *
         * Render the intersection.
         *
         */
        match intersection_result {
            IntersectResult::Intersect(points) => {
                draw_polyline(
                    points
                        .segments()
                        .map(|s| (mquad_from_na(s.a), mquad_from_na(s.b)))
                        .collect(),
                    Color::new(0f32, 1f32, 0f32, 1f32),
                );
                set_default_camera();
                easy_draw_text("Intersection found!");
            }
            IntersectResult::Negative => {
                set_default_camera();
                easy_draw_text("No intersection found, the shape is below the plane.");
            }
            IntersectResult::Positive => {
                set_default_camera();
                easy_draw_text("No intersection found, the shape is above the plane.");
            }
        }
        next_frame().await
    }
}
