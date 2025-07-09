//! Builds on top of `point_in_poly2d` example but uses [PointQuery] methods and different shapes.

mod common_macroquad2d;

use core::f32;

use common_macroquad2d::draw_point;
use macroquad::prelude::*;
use nalgebra::{Isometry, Point2, Rotation, Translation, UnitComplex};
use parry2d::query::gjk::GjkOptions;
use parry2d::shape::Compound;
use parry2d::shape::{ConvexPolygon, Shape, SharedShape};

const RENDER_SCALE: f32 = 30.0;

#[macroquad::main("points_in_poly2d")]
async fn main() {
    env_logger::init();
    let mut simple_convex = SharedShape::new(simple_convex());
    let mut simple_compound = Box::new(Compound::new(vec![(
        Isometry::default(),
        simple_convex.clone(),
    )]));
    let test_points = grid_points();

    let animation_rotation = UnitComplex::new(0.02);
    let polygon_render_pos = Point2::new(screen_width() / 2.0, screen_height() / 2.0);

    for i in 0.. {
        let shape_to_query = if (i) % 2 == 0 {
            simple_convex.0.as_ref()
        } else {
            (simple_compound.as_ref() as &dyn Shape)
        };

        clear_background(BLACK);

        // TODO: Display the shape
        /*         polygon
                    .iter_mut()
                    .for_each(|pt| *pt = animation_rotation * *pt);

                draw_polygon(&polygon, RENDER_SCALE, polygon_render_pos, BLUE);
        */
        /*
         * Compute polygon intersections.
         */
        for point in &test_points {
            let pos12 = Isometry::default().inv_mul(&Isometry::from_parts(
                Translation::from(point.coords),
                Rotation::identity(),
            ));
            if shape_to_query.contains_local_point(
                &(pos12 * point),
                // FIXME: Thierry (pr 298): options should contain everything necessary to be dispatched.
                &GjkOptions {
                    espilon_tolerance: f32::EPSILON * 10000.0,
                    nb_max_iterations: 100,
                },
            ) {
                draw_point(*point, RENDER_SCALE, polygon_render_pos, RED);
            } else {
                draw_point(*point, RENDER_SCALE, polygon_render_pos, GREEN);
            }
        }

        next_frame().await
    }
}

fn simple_convex() -> ConvexPolygon {
    let to_cast_against = ConvexPolygon::from_convex_polyline(
        [
            [-24.0, 0.0].into(),
            [0.0, -24.0].into(),
            [24.0, 0.0].into(),
            [0.0, 24.0].into(),
        ]
        .into(),
    )
    .unwrap();
    to_cast_against
}

fn grid_points() -> Vec<Point2<f32>> {
    let count = 80 * 5;
    let spacing = 0.6 / 5.0;
    let mut pts = vec![];
    for i in 0..count {
        for j in 0..count {
            pts.push(Point2::new(
                (i as f32 - count as f32 / 2.0) * spacing,
                (j as f32 - count as f32 / 2.0) * spacing,
            ));
        }
    }
    pts
}
