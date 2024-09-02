use macroquad::models::Vertex;
use macroquad::prelude::*;
use nalgebra::{Point3, UnitVector3, Vector3};
use parry3d::math::Real;
use parry3d::query::IntersectResult;
use parry3d::shape::{Cuboid, TriMesh};

#[macroquad::main("parry3d::query::PlaneIntersection")]
async fn main() {
    let trimesh = Cuboid::new(Vector3::repeat(1.0)).to_trimesh();

    let camera_pos = Vec3::new(-1.5f32, 2.5f32, -3f32);

    let mesh = mquad_mesh_from_points(&trimesh, camera_pos);
    let trimesh = TriMesh::new(trimesh.0, trimesh.1);

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
                draw_text("Intersection found!");
            }
            IntersectResult::Negative => {
                set_default_camera();
                draw_text("No intersection found, the shape is below the plane.");
            }
            IntersectResult::Positive => {
                set_default_camera();
                draw_text("No intersection found, the shape is above the plane.");
            }
        }
        next_frame().await
    }
}

fn mquad_mesh_from_points(trimesh: &(Vec<Point3<Real>>, Vec<[u32; 3]>), camera_pos: Vec3) -> Mesh {
    let (points, indices) = trimesh;
    // Transform the parry mesh into a mquad Mesh
    let (mquad_points, mquad_indices) = (
        points
            .iter()
            .map(|p| Vertex {
                position: mquad_from_na(*p),
                uv: Vec2::new(p.x, p.y),
                color: DARKGRAY.into(),
                normal: Vec4::ZERO,
            })
            .collect(),
        indices.iter().flatten().map(|v| *v as u16).collect(),
    );

    // Macroquad doesn´t support adding normals to vertices, so we'll bake a color into these vertices.
    // See https://github.com/not-fl3/macroquad/issues/321.

    // Compute the normal of each vertex, making them unique
    let vertices: Vec<Vertex> = mquad_compute_normals(&mquad_points, &mquad_indices, camera_pos);
    // Regenerate the index for each vertex.
    let indices: Vec<u16> = (0..vertices.len() * 3)
        .into_iter()
        .map(|i| i as u16)
        .collect();
    let mesh = Mesh {
        vertices,
        indices,
        texture: None,
    };
    mesh
}

fn mquad_compute_normals(points: &Vec<Vertex>, indices: &Vec<u16>, cam_pos: Vec3) -> Vec<Vertex> {
    let mut vertices: Vec<Vertex> = Vec::<Vertex>::new();
    for indices in indices.chunks(3) {
        let v0 = &points[indices[0] as usize];
        let v1 = &points[indices[1] as usize];
        let v2 = &points[indices[2] as usize];
        let normal = (v0.position - v2.position)
            .cross(v1.position - v2.position)
            .normalize();
        let brightness_mod = 0.2 + (0.8 / 2.) * (normal.dot(cam_pos) + 1.);

        for &i in indices.iter() {
            let mut color = points[i as usize].color;
            color[0] = (color[0] as f32 * brightness_mod) as u8;
            color[1] = (color[1] as f32 * brightness_mod) as u8;
            color[2] = (color[2] as f32 * brightness_mod) as u8;

            vertices.push(Vertex {
                position: points[i as usize].position,
                uv: Vec2::ZERO,
                color: color,
                normal: Vec4::ZERO,
            });
        }
    }
    vertices
}

fn draw_polyline(polygon: Vec<(Vec3, Vec3)>, color: Color) {
    for i in 0..polygon.len() {
        let a = polygon[i].0;
        let b = polygon[i].1;
        draw_line_3d(a, b, color);
    }
}

fn mquad_from_na(a: Point3<Real>) -> Vec3 {
    Vec3::new(a.x, a.y, a.z)
}

fn draw_text(text: &str) {
    macroquad::text::draw_text(text, 10.0, 48.0 + 18.0, 30.0, WHITE);
}
