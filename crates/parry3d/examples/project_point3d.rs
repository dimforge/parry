use macroquad::prelude::*;
use nalgebra::Vector3;
use parry3d::query::PointQuery;
use parry3d::shape::{Cuboid, TriMesh, TriMeshFlags};

mod common_macroquad3d;
use common_macroquad3d::*;

#[macroquad::main("project_point3d")]
async fn main() {
    let trimesh = Cuboid::new(Vector3::new(0.2, 0.5, 1.0)).to_trimesh();

    let mesh = mquad_mesh_from_points(
        &trimesh,
        Vec3::new(1f32, 3f32, 3f32),
        Color::from_rgba(200, 200, 200, 150),
    );
    let (points, indices) = trimesh;
    let trimesh = TriMesh::with_flags(points, indices, TriMeshFlags::ORIENTED).unwrap();
    for _i in 1.. {
        clear_background(BLACK);

        let elapsed_time = get_time() as f32;
        let slow_elapsed_time = elapsed_time / 3.0;

        let point_to_project = lissajous_3d(slow_elapsed_time);
        let projected_point = trimesh.project_local_point(&na_from_mquad(point_to_project), true);

        let slow_elapsed_time = slow_elapsed_time * 0.7;
        // Setup 3D camera.
        set_camera(&Camera3D {
            position: Vec3::new(
                slow_elapsed_time.sin() * 3.0,
                slow_elapsed_time.sin(),
                slow_elapsed_time.cos() * 3.0,
            ),
            up: Vec3::Y,
            target: Vec3::ZERO,
            ..Default::default()
        });

        /*
         *
         * Render the projection
         *
         */
        let color = if projected_point.is_inside {
            RED
        } else {
            YELLOW
        };

        draw_line_3d(
            point_to_project,
            mquad_from_na(projected_point.point),
            color,
        );
        draw_sphere(point_to_project, 0.1, None, color);

        draw_line_3d(
            point_to_project,
            mquad_from_na(projected_point.point),
            color,
        );

        // fixed point inside the shape
        let point_to_project = Vec3::ZERO;
        let projected_point = trimesh.project_local_point(&na_from_mquad(point_to_project), true);
        let color = if projected_point.is_inside {
            RED
        } else {
            YELLOW
        };
        draw_sphere(point_to_project, 0.1, None, color);

        draw_line_3d(
            point_to_project,
            mquad_from_na(projected_point.point),
            color,
        );
        // Mesh is rendered in the back, so we can see the other graphics elements
        draw_mesh(&mesh);

        next_frame().await
    }
}
