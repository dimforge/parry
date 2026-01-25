mod utils2d;

use kiss3d::prelude::*;
use parry2d::math::{Pose, Rotation};
use parry2d::query::PointQuery;
use parry2d::shape::{Cuboid, TriMesh, TriMeshFlags};
use utils2d::{draw_circle, draw_line_2d, draw_trimesh2, lissajous_2d};

#[kiss3d::main]
async fn main() {
    let mut window = Window::new("project_point2d").await;
    let mut camera = PanZoomCamera2d::new(Vec2::ZERO, 4.0);
    let mut scene = SceneNode2d::empty();

    let scale = 200f32;
    let (points, indices) = Cuboid::new(Vec2::new(0.2 * scale, 0.5 * scale)).to_trimesh();

    let trimesh = TriMesh::with_flags(points, indices, TriMeshFlags::ORIENTED).unwrap();

    let start_time = web_time::Instant::now();

    while window.render_2d(&mut scene, &mut camera).await {
        let elapsed_time = start_time.elapsed().as_secs_f32();
        let slow_elapsed_time = elapsed_time / 3.0;

        let offset = Vec2::ZERO;

        let point_to_project = lissajous_2d(slow_elapsed_time) * scale + offset;
        let projected_point = trimesh.project_point(
            &Pose::from_parts(offset, Rotation::identity()),
            point_to_project,
            true,
        );

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

        draw_line_2d(&mut window, point_to_project, projected_point.point, color);
        draw_circle(&mut window, point_to_project, 10f32, color);

        // fixed local point inside the shape
        let point_to_project = Vec2::ZERO;
        let projected_point = trimesh.project_local_point(point_to_project, true);
        let color = if projected_point.is_inside {
            RED
        } else {
            YELLOW
        };
        // convert to "world" space
        let point_to_project_world = point_to_project * scale + offset;
        draw_circle(&mut window, point_to_project_world, 10f32, color);

        draw_line_2d(
            &mut window,
            point_to_project_world,
            projected_point.point * scale + offset,
            color,
        );
        // Mesh is rendered in the back, so we can see the other graphics elements
        draw_trimesh2(&mut window, &trimesh, offset);
    }
}
