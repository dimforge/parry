use macroquad::models::Vertex;
use macroquad::prelude::*;
use nalgebra::{Point3, UnitVector3, Vector3};
use parry3d::math::{Isometry, Real};
use parry3d::query::{IntersectResult, PointQuery};
use parry3d::shape::TriMesh;

fn build_diamond(position: &Isometry<Real>) -> (Vec<Point3<Real>>, Vec<[u32; 3]>) {
    // Two tetrahedrons sharing a face
    let points = vec![
        position * Point3::new(0.0, 2.0, 0.0),
        position * Point3::new(-2.0, -1.0, 0.0),
        position * Point3::new(0.0, 0.0, 2.0),
        position * Point3::new(2.0, -1.0, 0.0),
        position * Point3::new(0.0, 0.0, -2.0),
    ];

    let indices = vec![
        [0u32, 1, 2],
        [0, 2, 3],
        [1, 2, 3],
        [0, 1, 4],
        [0, 4, 3],
        [1, 4, 3],
    ];

    (points, indices)
}

#[macroquad::main("parry3d::query::PlaneIntersection")]
async fn main() {
    //
    // This is useful to test for https://github.com/dimforge/parry/pull/248
    let _points = vec![
        Point3::from([0.0, 0.0, 0.0]),
        Point3::from([0.0, 0.0, 1.0]),
        Point3::from([1.0, 0.0, 0.0]),
        Point3::from([1.0, 0.0, 1.0]),
    ];
    let _indices: Vec<[u32; 3]> = vec![[0, 1, 2], [1, 3, 2]];
    //
    //

    let (points, indices) = build_diamond(&Isometry::identity());

    let mesh = Mesh {
        vertices: points
            .iter()
            .map(|p| Vertex {
                position: mquad_from_na(*p),
                uv: Vec2::new(p.x, p.y),
                color: Color::new(0.9, 0.9, 0.9, 0.7),
            })
            .collect(),
        indices: indices.iter().flatten().map(|v| *v as u16).collect(),
        texture: None,
    };
    let trimesh = TriMesh::new(points, indices);

    for _i in 1.. {
        clear_background(BLACK);

        let elapsed_time = get_time() as f32;
        let slow_elapsed_time = elapsed_time / 3.0;

        let sin = (elapsed_time / 3f32).sin();
        let bias = 1.5 * sin.abs();
        let rotation = Quat::from_axis_angle(
            Vec3::new(slow_elapsed_time.sin(), slow_elapsed_time.cos(), 1f32),
            (elapsed_time * 50f32).to_radians(),
        );
        let up_plane_vector = rotation * Vec3::Y;
        let point_to_project = up_plane_vector * bias;
        let projected_point = trimesh.project_point(
            &Isometry::identity(),
            &na_from_mquad(point_to_project),
            true,
        );

        let slow_elapsed_time = slow_elapsed_time / 2.0;
        // Going 3d!
        set_camera(&Camera3D {
            position: Vec3::new(
                slow_elapsed_time.sin() * 5.0,
                slow_elapsed_time.sin(),
                slow_elapsed_time.cos() * 5.0,
            ),
            up: Vec3::new(0f32, 1f32, 0f32),
            target: Vec3::new(0.5f32, 0f32, 0.5f32),
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

        // fixed point inside
        let point_to_project = Vec3::ZERO;
        let projected_point = trimesh.project_point(
            &Isometry::identity(),
            &na_from_mquad(point_to_project),
            true,
        );
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

        // Mesh is rendered in the back.
        draw_mesh(&mesh);

        next_frame().await
    }
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

fn na_from_mquad(a: Vec3) -> Point3<Real> {
    Point3::new(a.x, a.y, a.z)
}
