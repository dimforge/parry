mod common_macroquad2d;

use common_macroquad2d::{draw_line_2d, draw_trimesh2, lissajous_2d, mquad_from_na, na_from_mquad};
use macroquad::prelude::*;
use nalgebra::{Point3, UnitComplex, Vector2};
use parry2d::math::{Isometry, Translation};
use parry2d::query::PointQuery;
use parry2d::shape::{Cuboid, TriMesh, TriMeshFlags};

#[macroquad::main("project_point2d")]
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

    let trimesh = TriMesh::with_flags(points, indices, TriMeshFlags::ORIENTED).unwrap();
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
