mod utils3d;

use core::f32::consts::{FRAC_PI_2, FRAC_PI_4, FRAC_PI_6};

use kiss3d::prelude::*;
use parry3d::transformation;
use utils3d::{create_mesh_from_trimesh, lissajous_3d_with_params};

#[kiss3d::main]
async fn main() {
    let mut window = Window::new("convex_hull3d").await;
    let mut camera = OrbitCamera3d::new(Vec3::new(8.0, 8.0, 8.0), Vec3::new(0.5, 0.0, 0.5));
    let mut scene = SceneNode3d::empty();
    scene
        .add_light(Light::point(100.0))
        .set_position(Vec3::new(5.0, 10.0, 3.0));

    let count = 9;
    let mut pts = vec![Vec3::ZERO; count];

    let start_time = web_time::Instant::now();
    let mut mesh_node: Option<SceneNode3d> = None;

    while window.render_3d(&mut scene, &mut camera).await {
        let elapsed_time = start_time.elapsed().as_secs_f32();
        let elapsed_time_slow = elapsed_time * 0.2;

        for (i, pt) in pts.iter_mut().enumerate() {
            *pt = lissajous_3d_with_params(
                (i * i) as f32 + elapsed_time_slow,
                2.0 + i as f32 / 3.0,
                1f32 + (i as f32).sin() * 0.2,
                (i as f32 / count as f32) + elapsed_time_slow.cos() * 0.1,
                (elapsed_time_slow as f32 + i as f32).cos() * 0.1 + FRAC_PI_2,
                FRAC_PI_4,
                FRAC_PI_6,
            ) * 5f32;
            window.draw_point(*pt, RED, 20.0);
        }

        /*
         *
         * Compute the convex hull.
         *
         */
        let convex_hull = transformation::convex_hull(&pts);

        // Remove previous mesh if exists
        if let Some(mut node) = mesh_node.take() {
            node.remove();
        }

        // Create new mesh from convex hull
        let mesh = create_mesh_from_trimesh(convex_hull.0, convex_hull.1);
        mesh_node = Some(
            scene
                .add_mesh(mesh, Vec3::splat(1.0))
                .set_color(LIGHT_GREEN),
        );
    }
}
