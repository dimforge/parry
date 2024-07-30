use macroquad::models::Vertex;
use macroquad::prelude::*;
use nalgebra::{Point3, UnitVector3, Vector3};
use parry3d::math::{Isometry, Real};
use parry3d::query::IntersectResult;
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
                position: Vec3::new(p.x, p.y, p.z),
                uv: Vec2::new(p.x, p.y),
                color: WHITE,
            })
            .collect(),
        indices: indices.as_flattened().iter().map(|v| *v as u16).collect(),
        texture: None,
    };
    let trimesh = TriMesh::new(points, indices);

    for _i in 1.. {
        clear_background(BLACK);

        let elapsed_time = get_time();

        let bias = -2.0 * (elapsed_time as f32 / 3f32).sin();
        let rotation = Quat::from_axis_angle(Vec3::Z, (elapsed_time as f32 * 50f32).to_radians());
        let up_plane_vector = rotation * Vec3::Y;
        let intersection_result = trimesh.intersection_with_local_plane(
            &UnitVector3::new_normalize(Vector3::<Real>::new(
                up_plane_vector.x,
                up_plane_vector.y,
                up_plane_vector.z,
            )),
            bias,
            0.0005,
        );

        // Going 3d!
        set_camera(&Camera3D {
            position: Vec3::new(0f32, 3f32, -3f32),
            up: Vec3::new(0f32, 1f32, 0f32),
            target: Vec3::new(0.5f32, 0f32, 0.5f32),
            ..Default::default()
        });

        let plane_center = up_plane_vector * bias;
        draw_line_3d(plane_center, plane_center + up_plane_vector, GREEN);
        draw_mesh(&mesh);
        draw_grid_ex(10, 0.333, BLUE, RED, plane_center, rotation);

        /*
         *
         * Render the intersection
         *
         */
        match intersection_result {
            IntersectResult::Intersect(points) => {
                draw_polyline(
                    points
                        .segments()
                        .map(|s| {
                            (
                                Vec3::new(s.a.x, s.a.y, s.a.z),
                                Vec3::new(s.b.x, s.b.y, s.b.z),
                            )
                        })
                        .collect(),
                    Color::new(0f32, 1f32, 0f32, 1f32),
                );
            }
            IntersectResult::Negative => {
                set_default_camera();
                draw_text(
                    format!("No intersection found, the plane goes below the shape.").as_str(),
                    10.0,
                    48.0 + 18.0,
                    30.0,
                    WHITE,
                );
            }
            IntersectResult::Positive => {
                set_default_camera();
                draw_text(
                    format!("No intersection found, the plane goes above the shape.").as_str(),
                    10.0,
                    48.0 + 18.0,
                    30.0,
                    WHITE,
                );
            }
        }
        next_frame().await
    }
}

fn draw_polyline(polygon: Vec<(Vec3, Vec3)>, color: Color) {
    for i in 0..polygon.len() {
        let a = polygon[i].0;
        let b = polygon[i].1;
        draw_line_3d(Vec3::new(a.x, a.y, a.z), Vec3::new(b.x, b.y, b.z), color);
    }
}
