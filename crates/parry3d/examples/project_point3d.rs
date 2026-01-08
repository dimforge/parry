mod utils3d;

use kiss3d::prelude::*;
use parry3d::query::PointQuery;
use parry3d::shape::{Cuboid, TriMesh, TriMeshFlags};
use utils3d::{create_mesh_from_trimesh, lissajous_3d};

#[kiss3d::main]
async fn main() {
    let mut window = Window::new("project_point3d").await;
    let mut scene = SceneNode3d::empty();
    scene
        .add_light(Light::point(100.0))
        .set_position(Vec3::new(1.0, 3.0, 3.0));

    let (points, indices) = Cuboid::new(Vec3::new(0.2, 0.5, 1.0)).to_trimesh();
    let trimesh =
        TriMesh::with_flags(points.clone(), indices.clone(), TriMeshFlags::ORIENTED).unwrap();

    // Add the mesh to the scene
    let mesh = create_mesh_from_trimesh(points, indices);
    scene
        .add_mesh(mesh, Vec3::splat(1.0))
        .set_color(Color::new(0.8, 0.8, 0.8, 0.6));

    let start_time = web_time::Instant::now();
    let mut camera = OrbitCamera3d::new(Vec3::new(3.0, 1.0, 3.0), Vec3::ZERO);

    while window.render_3d(&mut scene, &mut camera).await {
        let elapsed_time = start_time.elapsed().as_secs_f32();
        let slow_elapsed_time = elapsed_time / 3.0;

        let point_to_project = lissajous_3d(slow_elapsed_time);
        let projected_point = trimesh.project_local_point(point_to_project, true);

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

        window.draw_line(point_to_project, projected_point.point, color, 2.0, false);
        window.draw_point(point_to_project, color, 10.0);

        // fixed point inside the shape
        let point_to_project = Vec3::ZERO;
        let projected_point = trimesh.project_local_point(point_to_project, true);
        let color = if projected_point.is_inside {
            RED
        } else {
            YELLOW
        };
        window.draw_point(point_to_project, color, 10.0);
        window.draw_line(point_to_project, projected_point.point, color, 2.0, false);
    }
}
