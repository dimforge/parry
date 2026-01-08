mod utils3d;

use kiss3d::prelude::*;
use parry3d::query::IntersectResult;
use parry3d::shape::{Cuboid, TriMesh};
use utils3d::{create_mesh_from_trimesh, draw_polyline, draw_text};

#[kiss3d::main]
async fn main() {
    let mut window = Window::new("plane_intersection").await;
    let mut camera = OrbitCamera3d::new(Vec3::new(-1.5, 2.5, -3.0), Vec3::new(0.5, 0.0, 0.5));
    let mut scene = SceneNode3d::empty();
    scene
        .add_light(Light::point(100.0))
        .set_position(Vec3::new(-1.0, 3.5, -3.0));

    let font = Font::default();

    let (points, indices) = Cuboid::new(Vec3::splat(1.0)).to_trimesh();
    let trimesh = TriMesh::new(points.clone(), indices.clone()).unwrap();

    // Add mesh to scene
    let mesh = create_mesh_from_trimesh(points, indices);
    scene
        .add_mesh(mesh, Vec3::splat(1.0))
        .set_color(Color::new(0.3, 0.3, 0.3, 1.0));

    let start_time = web_time::Instant::now();

    while window.render_3d(&mut scene, &mut camera).await {
        let elapsed_time = start_time.elapsed().as_secs_f64();

        // Animated rotation for the intersection plane.
        let bias = -1.2 * (elapsed_time as f32 / 3f32).sin();
        let rotation = Quat::from_axis_angle(Vec3::Z, (elapsed_time as f32 * 40f32).to_radians());
        let up_plane_vector = rotation * Vec3::Y;

        // Get the intersection polyline.
        let intersection_result =
            trimesh.intersection_with_local_plane(up_plane_vector.normalize(), bias, 0.0005);

        // Draw the plane normal
        let plane_center = up_plane_vector * bias;
        window.draw_line(
            plane_center,
            plane_center + up_plane_vector,
            GREEN,
            4.0,
            false,
        );

        // Draw a grid to represent the plane
        draw_grid(&mut window, plane_center, rotation, 10, 0.333);

        /*
         *
         * Render the intersection.
         *
         */
        match intersection_result {
            IntersectResult::Intersect(points) => {
                draw_polyline(
                    &mut window,
                    points.segments().map(|s| (s.a, s.b)).collect(),
                    Color::new(0f32, 1f32, 0f32, 1f32),
                );
                draw_text(&mut window, &font, "Intersection found!");
            }
            IntersectResult::Negative => {
                draw_text(
                    &mut window,
                    &font,
                    "No intersection found, the shape is below the plane.",
                );
            }
            IntersectResult::Positive => {
                draw_text(
                    &mut window,
                    &font,
                    "No intersection found, the shape is above the plane.",
                );
            }
        }
    }
}

fn draw_grid(window: &mut Window, center: Vec3, rotation: Quat, subdivisions: i32, spacing: f32) {
    let half_size = subdivisions as f32 * spacing / 2.0;

    for i in 0..=subdivisions {
        let offset = -half_size + i as f32 * spacing;

        // Lines along X axis
        let start = rotation * Vec3::new(-half_size, 0.0, offset) + center;
        let end = rotation * Vec3::new(half_size, 0.0, offset) + center;
        let color = if i == subdivisions / 2 { RED } else { BLUE };
        window.draw_line(start, end, color, 2.0, false);

        // Lines along Z axis
        let start = rotation * Vec3::new(offset, 0.0, -half_size) + center;
        let end = rotation * Vec3::new(offset, 0.0, half_size) + center;
        let color = if i == subdivisions / 2 { RED } else { BLUE };
        window.draw_line(start, end, color, 2.0, false);
    }
}
