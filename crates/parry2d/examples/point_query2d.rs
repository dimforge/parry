//! Builds on top of `point_in_poly2d` example but uses [PointQuery] methods and different shapes.

mod common_macroquad2d;
use core::f32;

use common_macroquad2d::draw_point;
use macroquad::prelude::*;
use nalgebra::{Isometry, Point2, Rotation, Translation};
use parry2d::query::details::point_query::QueryOptions;
use parry2d::query::gjk::GjkOptions;
use parry2d::query::point::QueryOptionsDispatcherMap;
use parry2d::shape::Compound;
use parry2d::shape::{ConvexPolygon, Shape, SharedShape};

use crate::common_macroquad2d::easy_draw_text;

const RENDER_SCALE: f32 = 30.0;

#[macroquad::main("points_in_poly2d")]
async fn main() {
    env_logger::init();
    let simple_convex = SharedShape::new(simple_convex());
    let simple_compound = Box::new(Compound::new(vec![(
        Isometry::default(),
        simple_convex.clone(),
    )]));
    let test_points = grid_points();

    let polygon_render_pos = Point2::new(screen_width() / 2.0, screen_height() / 2.0);

    // Create an initial dispatcher
    let mut query_options_dispatcher = QueryOptionsDispatcherMap::default();

    for i in 0.. {
        clear_background(BLACK);

        let gjk_options = query_options_dispatcher
            .get_option_mut::<GjkOptions>()
            .unwrap();

        // loops from default epsilon to an arbitrarily chosen slightly higher value.
        gjk_options.epsilon_tolerance =
            f32::EPSILON + (((i as f32 / 10f32).sin() + 1f32) / 2f32) * 0.000002f32;
        let (shape_to_query, options) = if (i) % 2 == 0 {
            (simple_convex.0.as_ref(), &*gjk_options as &dyn QueryOptions)
        } else {
            (
                simple_compound.as_ref() as &dyn Shape,
                &query_options_dispatcher as &dyn QueryOptions,
            )
        };
        for point in &test_points {
            let pos12 = Isometry::default().inv_mul(&Isometry::from_parts(
                Translation::from(point.coords),
                Rotation::identity(),
            ));
            if shape_to_query.contains_local_point(&(pos12 * point), options) {
                draw_point(*point, RENDER_SCALE, polygon_render_pos, GREEN);
            } else {
                draw_point(*point, RENDER_SCALE, polygon_render_pos, DARKGRAY);
            }
        }
        let gjk_options = query_options_dispatcher.get_option::<GjkOptions>().unwrap();

        easy_draw_text(&format!("tolerance: {:.7}", gjk_options.epsilon_tolerance));

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
    let count = 50;
    let spacing = 0.4;
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
