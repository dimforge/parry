use std::f32::consts::{FRAC_PI_2, FRAC_PI_4, FRAC_PI_6};

use macroquad::gizmos::{draw_gizmos, gizmos_add_line, init_gizmos};
use macroquad::math::{vec2, vec3, Mat4, Vec2, Vec3};
use macroquad::quad_gl::camera::{Camera, Projection};
use macroquad::quad_gl::models::CpuMesh;
use macroquad::quad_gl::QuadGl;
use macroquad::quad_gl::{
    camera::Environment,
    color::{self},
    scene::Shader,
};
use macroquad::window::{next_frame, screen_height, screen_width};

use nalgebra::{Point3, Vector3};
use parry3d::math::{Isometry, Real};
use parry3d::query::PointQuery;
use parry3d::shape::{Cuboid, TriMesh, TriMeshFlags};

fn lissajous_3d(t: f32) -> Vec3 {
    // Some hardcoded parameters to have a pleasing lissajous trajectory.
    let (a, b, c, delta_x, delta_y, delta_z) = (1.0, 3.0, 2.0, FRAC_PI_4, FRAC_PI_2, FRAC_PI_4);

    let x = (a * t + delta_x).sin();
    let y = (b * t + delta_y).sin();
    let z = (c * t + delta_z).sin();
    Vec3::new(x, y, z) * 0.8f32
}

fn main() {
    macroquad::start(Default::default(), |ctx| main_loop(ctx));
}

async fn main_loop(ctx: macroquad::Context) {
    init_gizmos(&ctx);

    let (points, indices) = Cuboid::new(Vector3::new(0.5, 0.5, 0.8)).to_trimesh();

    let mut scene = ctx.new_scene();

    let quad_gl = QuadGl::new(ctx.quad_ctx.clone());

    let cpu_mesh = mquad_mesh_from_parry(&indices, &points);

    let mut mesh = quad_gl.mesh(cpu_mesh, None);

    mesh.nodes[0].materials[0].shader = Shader::new(
        ctx.quad_ctx.lock().unwrap().as_mut(),
        vec![],
        Some(FRAGMENT),
        None,
    );

    scene.add_model(&mesh);
    let trimesh = TriMesh::with_flags(points, indices, TriMeshFlags::ORIENTED);
    let mut canvas = ctx.new_canvas();

    let mut camera = Camera {
        environment: Environment::SolidColor(color::BLACK),
        depth_enabled: true,
        projection: Projection::Perspective,
        position: vec3(0., 1.5, 4.),
        up: vec3(0., 1., 0.),
        target: vec3(0., 0., 0.),
        z_near: 0.1,
        z_far: 1500.0,
        ..Default::default()
    };

    // FIXME: that's framerate dependent
    let mut elapsed_time = 0.0f32;
    for _i in 1.. {
        elapsed_time += 0.1;
        ctx.clear_screen(color::BLACK);
        canvas.clear();

        let slow_elapsed_time = elapsed_time / 10.0;

        let point_to_project = lissajous_3d(slow_elapsed_time);
        let projected_point = trimesh.project_point(
            &Isometry::identity(),
            &na_from_mquad(point_to_project),
            true,
        );

        camera.position = Vec3::new(
            slow_elapsed_time.sin() * 5.0,
            slow_elapsed_time.sin() * 1.5,
            slow_elapsed_time.cos() * 5.0,
        );

        //camera.position = Vec3::new(-1.0, 2.5, -4.0);

        /*
         *
         * Render the projection
         *
         */

        let color = if projected_point.is_inside {
            color::RED
        } else {
            color::YELLOW
        };

        gizmos_add_line(
            false,
            point_to_project,
            mquad_from_na(projected_point.point),
        );
        // not working
        //let point = camera
        //    .proj_view()
        //    .0
        //    .project_point3(point_to_project.xy);
        //let point = vec2(
        //    (point.x / 2. + 0.5) * screen_width(),
        //    (0.5 - point.y / 2.) * screen_height(),
        //);

        //canvas.draw_circle(point.x, point.y, 10.0, color);
        //draw_sphere(point_to_project, 0.1, None, color);

        gizmos_add_line(
            false,
            point_to_project,
            mquad_from_na(projected_point.point),
        );

        // fixed point inside
        let point_to_project = Vec3::ZERO;
        let projected_point = trimesh.project_point(
            &Isometry::identity(),
            &na_from_mquad(point_to_project),
            true,
        );
        let color = if projected_point.is_inside {
            color::RED
        } else {
            color::YELLOW
        };
        //draw_sphere(point_to_project, 0.1, None, color);

        gizmos_add_line(
            false,
            point_to_project,
            mquad_from_na(projected_point.point),
        );

        ctx.root_ui().draw(&mut canvas);

        scene.draw(&camera);
        draw_gizmos(&camera);

        canvas.draw();

        next_frame().await
    }
}

fn mquad_mesh_from_parry(indices: &Vec<[u32; 3]>, points: &Vec<Point3<Real>>) -> CpuMesh {
    let m_indices = indices.iter().flatten().map(|v| *v as u16).collect();
    let m_vertices = points.iter().map(|p| mquad_from_na(*p)).collect();

    let mesh = compute_mesh_with_normals_per_face(&m_vertices, &m_indices);

    let nb_vertices = mesh.vertices.len();
    let cpu_mesh = CpuMesh(
        mesh.vertices,
        vec![vec2(0.0, 0.0); nb_vertices],
        mesh.normals,
        mesh.indices,
    );
    cpu_mesh
}

fn mquad_from_na(a: Point3<Real>) -> Vec3 {
    Vec3::new(a.x, a.y, a.z)
}

fn na_from_mquad(a: Vec3) -> Point3<Real> {
    Point3::new(a.x, a.y, a.z)
}

pub struct MeshData {
    pub vertices: Vec<Vec3>,
    pub indices: Vec<u16>,
    pub normals: Vec<Vec3>,
}

pub fn compute_mesh_with_normals_per_face(vertices: &Vec<Vec3>, indices: &Vec<u16>) -> MeshData {
    let mut result_vertices: Vec<Vec3> = Vec::<Vec3>::new();
    let mut normals: Vec<Vec3> = Vec::<Vec3>::new();
    for indices in indices.chunks(3) {
        let v0 = vertices[indices[0] as usize];
        let v1 = vertices[indices[1] as usize];
        let v2 = vertices[indices[2] as usize];

        let normal = (v0 - v2).cross(v1 - v2).normalize();

        for &i in indices.iter() {
            result_vertices.push(vertices[i as usize]);
            normals.push(normal);
        }
    }
    let vertices_count = result_vertices.len();
    MeshData {
        vertices: result_vertices,
        indices: (0..vertices_count).into_iter().map(|i| i as u16).collect(),
        normals,
    }
}

const VERTEX: &str = r#"
#include "common_vertex.glsl"

void vertex() {
}
"#;

const FRAGMENT: &str = r#"
varying vec3 out_normal;

void main() {
    vec3 norm = normalize(out_normal);
    vec3 lightDir = normalize(vec3(1.0, -1.0, 0.5));
    float diff = max(dot(norm, lightDir), 0.0);
    vec3 diffuse = diff * vec3(1.0) + vec3(0.3);

    gl_FragColor = vec4(diffuse,0.5);
}
"#;
