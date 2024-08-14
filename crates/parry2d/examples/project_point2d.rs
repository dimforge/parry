use std::f32::consts::{FRAC_PI_2, FRAC_PI_4, FRAC_PI_6};

use macroquad::models::Vertex;
use macroquad::prelude::*;
use nalgebra::{Point2, Point3, UnitComplex, UnitQuaternion, Vector2};
use parry2d::math::{Isometry, Real, Translation};
use parry2d::query::PointQuery;
use parry2d::shape::{Cuboid, TriMesh, TriMeshFlags};

fn lissajous_2d(t: f32) -> Vec2 {
    // Some hardcoded parameters to have a pleasing lissajous trajectory.
    let (a, b, delta_x, delta_y) = (3.0, 2.0, FRAC_PI_2, FRAC_PI_4);

    let x = (a * t + delta_x).sin();
    let y = (b * t + delta_y).sin();
    Vec2::new(x, y) * 0.75f32
}

const VIRTUAL_WIDTH: f32 = 1280.0;
const VIRTUAL_HEIGHT: f32 = 720.0;

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

    let scale = 200f32;
    let (points, indices) = Cuboid::new(Vector2::new(0.2 * scale, 0.5 * scale)).to_trimesh();

    let trimesh = TriMesh::with_flags(points, indices, TriMeshFlags::ORIENTED);
    for _i in 1.. {
        clear_background(BLACK);

        let elapsed_time = get_time() as f32;
        let slow_elapsed_time = elapsed_time / 3.0;

        let offset = Vec2::new(screen_width() / 2f32, screen_height() / 2f32);

        let point_to_project = lissajous_2d(slow_elapsed_time) * scale + offset;
        let translation = Translation::new(offset.x, offset.y);
        let rot = UnitComplex::identity();
        let projected_point = trimesh.project_point(
            &Isometry::from_parts(translation, rot),
            &na_from_mquad(point_to_project),
            true,
        );

        let slow_elapsed_time = slow_elapsed_time * 0.7;

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

        draw_line_2d(
            point_to_project,
            mquad_from_na(projected_point.point),
            color,
        );
        draw_circle(point_to_project.x, point_to_project.y, 10f32, color);

        draw_line_2d(
            point_to_project,
            mquad_from_na(projected_point.point),
            color,
        );

        // fixed local point inside the shape
        let point_to_project = Vec2::ZERO;
        let projected_point = trimesh.project_local_point(&na_from_mquad(point_to_project), true);
        let color = if projected_point.is_inside {
            RED
        } else {
            YELLOW
        };
        // convert to "world" space
        let point_to_project = point_to_project * scale + offset;
        draw_circle(point_to_project.x, point_to_project.y, 10f32, color);

        draw_line_2d(
            point_to_project,
            mquad_from_na(projected_point.point) * scale + offset,
            color,
        );
        // Mesh is rendered in the back, so we can see the other graphics elements
        draw_trimesh2(&trimesh, offset);

        next_frame().await
    }
}

fn mquad_from_na(a: Point2<Real>) -> Vec2 {
    Vec2::new(a.x, a.y)
}

fn na_from_mquad(a: Vec2) -> Point2<Real> {
    Point2::new(a.x, a.y)
}

fn draw_line_2d(a: Vec2, b: Vec2, color: Color) {
    draw_line(a.x, a.y, b.x, b.y, 2f32, color);
}

fn draw_trimesh2(trimesh: &TriMesh, offset: Vec2) {
    let vertices = trimesh.vertices();
    for v in trimesh.indices() {
        let v0 = mquad_from_na(vertices[v[0] as usize]) + offset;
        let v1 = mquad_from_na(vertices[v[1] as usize]) + offset;
        let v2 = mquad_from_na(vertices[v[2] as usize]) + offset;

        draw_line(v0.x, v0.y, v1.x, v1.y, 2f32, WHITE);
        draw_line(v0.x, v0.y, v2.x, v2.y, 2f32, WHITE);
        draw_line(v2.x, v2.y, v1.x, v1.y, 2f32, WHITE);
    }
}
