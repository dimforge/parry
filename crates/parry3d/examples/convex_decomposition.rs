mod utils3d;

use kiss3d::prelude::*;
use parry3d::{
    math::Real,
    shape::{SharedShape, TriMesh, TriMeshFlags},
};
use utils3d::{create_mesh_from_trimesh, draw_text, hue_to_rgb};

#[kiss3d::main]
async fn main() {
    let mut window = Window::new("convex_decomposition").await;
    let mut scene = SceneNode3d::empty();
    scene
        .add_light(Light::point(100.0))
        .set_position(Vec3::new(1.0, 3.0, 3.0));

    let font = Font::default();

    /*
     * Initialize the shapes.
     */
    let ::obj::Obj {
        data: ::obj::ObjData {
            position, objects, ..
        },
        ..
    } = ::obj::Obj::load("assets/tests/low_poly_bunny.obj").unwrap();

    let bunny_mesh = TriMesh::with_flags(
        position
            .iter()
            .map(|v| Vec3::new(v[0] as Real, v[1] as Real, v[2] as Real))
            .collect::<Vec<_>>(),
        objects[0].groups[0]
            .polys
            .iter()
            .map(|p| [p.0[0].0 as u32, p.0[1].0 as u32, p.0[2].0 as u32])
            .collect::<Vec<_>>(),
        TriMeshFlags::all(),
    )
    .unwrap();

    // Show loading message
    draw_text(
        &mut window,
        &font,
        "Please wait while convex decomposition is being computed...",
    );

    // Dummy render to display the message
    {
        let mut camera = OrbitCamera3d::new(Vec3::new(5.5, 3.0, 5.5), Vec3::new(0.0, 1.0, 0.0));
        window.render_3d(&mut scene, &mut camera).await;
    }

    let mesh_vertices = bunny_mesh.vertices();
    let mesh_indices = bunny_mesh.indices();
    let convex_mesh = SharedShape::convex_decomposition(&mesh_vertices, &mesh_indices);
    let trimesh_convex_compound = convex_mesh.as_compound().unwrap();

    let shapes_count = trimesh_convex_compound.shapes().len() as u32;

    // Add all convex parts to the scene
    for (i, s) in trimesh_convex_compound.shapes().iter().enumerate() {
        let (vertices, indices) = s.1.as_convex_polyhedron().unwrap().to_trimesh();
        let mesh = create_mesh_from_trimesh(vertices, indices);

        let (r, g, b) = hue_to_rgb(i as f32 / 6 as f32);
        scene
            .add_mesh(mesh, Vec3::splat(1.0))
            .set_color(Color::new(r, g, b, 1.0));
    }

    let mut camera = OrbitCamera3d::new(Vec3::new(5.5, 3.0, 5.5), Vec3::new(0.0, 1.0, 0.0));

    while window.render_3d(&mut scene, &mut camera).await {
        draw_text(
            &mut window,
            &font,
            &format!("Number of shapes: {}", shapes_count),
        );
    }
}
