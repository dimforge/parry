use nalgebra::{Point2, UnitComplex};
use parry2d::{self as parry, utils::point_in_poly2d};

use bevy::{color::palettes, prelude::*};

type Point = Point2<f32>;

#[derive(Resource)]
pub struct RapierContext {
    pub spikes: Vec<Point>,
    pub test_points: Vec<Vec2>,
}

fn main() {
    let mut app = App::new();

    app.add_plugins(DefaultPlugins);
    app.add_systems(Startup, init_physics);
    app.add_systems(Update, (rotate, show_physics).chain());

    app.run();
}

fn init_physics(mut commands: Commands) {
    let spikes = spikes_polygon();
    let test_points = grid_points();

    commands.insert_resource(RapierContext {
        spikes,
        test_points,
    });
    commands.spawn(Camera2dBundle::default());
}

fn rotate(mut physics: ResMut<RapierContext>) {
    let animation_rotation = UnitComplex::new(0.02);
    physics
        .spikes
        .iter_mut()
        .for_each(|pt| *pt = animation_rotation * *pt);
}

fn show_physics(mut gizmos: Gizmos, physics: Res<RapierContext>) {
    for p in &physics.test_points {
        let color = if point_in_poly2d(&Point::from(*p), &physics.spikes) {
            palettes::basic::RED
        } else {
            palettes::basic::GREEN
        };
        gizmos.circle_2d(*p, 3f32, color);
    }
    let spike_len = physics.spikes.len();
    for i in 0..spike_len {
        let a = physics.spikes[i];
        let b = physics.spikes[(i + 1) % spike_len];
        gizmos.line_2d(a.into(), b.into(), palettes::basic::BLUE);
    }
}

fn spikes_polygon() -> Vec<Point> {
    let teeths = 3;
    let width = 15.0 * 30f32;
    let height = 7.5 * 30f32;
    let tooth_width = width / (teeths as f32);
    let center = Point::new(width / 2.0, height / 2.0);

    let mut polygon: Vec<Point> = vec![
        (Point::new(width, 0.0) - center).into(),
        (Point::new(width, height) - center).into(),
        (Point::new(0.0, height) - center).into(),
    ];

    for i in 0..teeths {
        let x = i as f32 * tooth_width;
        polygon.push((Point::new(x, 0.0) - center).into());
        polygon.push((Point::new(x + tooth_width / 2.0, height * 1.5) - center).into());
    }

    polygon
}

fn grid_points() -> Vec<Vec2> {
    let count = 40;
    let spacing = 20f32;
    let mut pts = vec![];
    for i in 0..count {
        for j in 0..count {
            pts.push(Vec2::new(
                (i as f32 - count as f32 / 2.0) * spacing,
                (j as f32 - count as f32 / 2.0) * spacing,
            ));
        }
    }
    pts
}
