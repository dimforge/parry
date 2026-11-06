// Shape-casting against a `HeightField` skips the triangles the moving shape can't hit within
// `max_time_of_impact`, or before the nearest hit found so far. These check that the result is
// still the nearest of the hits of a shape-cast against each triangle on its own.

use approx::assert_relative_eq;
use parry3d::math::{Pose, Real, Rotation, Vector};
use parry3d::query::details::cast_shapes_heightfield_shape;
use parry3d::query::{
    self, ClosestPoints, Contact, DefaultQueryDispatcher, NonlinearRigidMotion, QueryDispatcher,
    ShapeCastHit, ShapeCastOptions, ShapeDistance, ShapeIntersection, Unsupported,
};
use parry3d::shape::{
    Ball, Capsule, Cone, Cuboid, HeightField, HeightFieldFlags, RoundCuboid, Shape, Triangle,
};
use parry3d::utils::Array2;
use std::sync::Mutex;

/// A heightfield `size` wide on a grid of `samples` × `samples`, with the given height function.
fn heightfield(samples: usize, size: Real, height: impl Fn(Real, Real) -> Real) -> HeightField {
    let heights = (0..samples)
        .flat_map(|i| (0..samples).map(move |j| (i as Real, j as Real)))
        .map(|(i, j)| height(i, j))
        .collect();
    HeightField::with_flags(
        Array2::new(samples, samples, heights),
        Vector::new(size, 1.0, size),
        HeightFieldFlags::FIX_INTERNAL_EDGES,
    )
}

/// The shape-casts against each triangle of `heightfield` on its own that hit.
fn triangle_hits(
    pos1: &Pose,
    heightfield: &HeightField,
    pos2: &Pose,
    vel2: Vector,
    g2: &dyn Shape,
    options: ShapeCastOptions,
) -> Vec<ShapeCastHit> {
    heightfield
        .triangles()
        .filter_map(|tri| {
            query::cast_shapes(pos1, Vector::ZERO, &tri, pos2, vel2, g2, options).unwrap()
        })
        .collect()
}

/// Shape-casts `g2` from `pos2` along `vel2` against `heightfield` at `pos1`, asserts that the
/// result is the nearest hit of the shape-casts against each triangle on its own (any one among
/// exact ties), and returns whether it hits.
fn assert_nearest_triangle_hit(
    pos1: &Pose,
    heightfield: &HeightField,
    pos2: &Pose,
    vel2: Vector,
    g2: &dyn Shape,
    options: ShapeCastOptions,
) -> bool {
    let hit = query::cast_shapes(pos1, Vector::ZERO, heightfield, pos2, vel2, g2, options).unwrap();
    let hits = triangle_hits(pos1, heightfield, pos2, vel2, g2, options);
    let nearest = hits.iter().map(|hit| hit.time_of_impact).reduce(Real::min);
    let context = || format!("{pos2:?}, {vel2:?}, {options:?}");

    assert_eq!(hit.map(|hit| hit.time_of_impact), nearest, "{}", context());
    if let Some(hit) = hit {
        // `ShapeCastHit` has no `PartialEq`; `Debug` prints each float so that it reads back
        // exactly.
        let hit = format!("{hit:?}");
        assert!(
            hits.iter().any(|tri_hit| format!("{tri_hit:?}") == hit),
            "{hit} is no triangle's hit: {}",
            context()
        );
    }
    hit.is_some()
}

fn shapes() -> [Box<dyn Shape>; 5] {
    [
        Box::new(Ball::new(0.4)),
        Box::new(Cuboid::new(Vector::new(0.6, 0.3, 0.4))),
        Box::new(Capsule::new_y(0.5, 0.4)),
        Box::new(Cone::new(0.4, 0.3)),
        Box::new(RoundCuboid {
            inner_shape: Cuboid::new(Vector::new(0.3, 0.2, 0.25)),
            border_radius: 0.15,
        }),
    ]
}

fn options(
    max_time_of_impact: Real,
    target_distance: Real,
    stop_at_penetration: bool,
) -> ShapeCastOptions {
    ShapeCastOptions {
        max_time_of_impact,
        target_distance,
        stop_at_penetration,
        compute_impact_geometry_on_penetration: true,
    }
}

/// Shape-casts over rolling and rough ground, standing on it, hovering over it and sunk in it,
/// down onto it, along it and up from it, short and long, with and without a target distance.
#[test]
fn heightfield_shape_cast_gives_the_nearest_triangle_hit() {
    // Waves about 10 m long and 0.5 m high, and ground that rises and falls by up to about 1 m
    // from one sample to the next.
    let grounds = [
        heightfield(11, 10.0, |i, j| 0.25 * ((0.6 * i).cos() + (0.6 * j).cos())),
        heightfield(11, 10.0, |i, j| 0.6 * ((1.3 * i).sin() + (1.1 * j).cos())),
    ];
    // The heightfield is posed off the origin and turned, so that the cast runs in its frame.
    let pos1 = Pose::from_parts(Vector::new(3.0, -1.0, 2.0), Rotation::from_rotation_y(0.3));
    // In the heightfield's frame. The traversal of the heightfield misses triangles on its own,
    // before any of them is skipped, when the shape moves along +X or +Z, or with no vertical
    // motion, so these don't.
    let directions = [
        Vector::new(-1.0, -0.01, 0.0),
        Vector::new(-1.0, 0.02, -1.0),
        Vector::new(-0.1, -0.03, -1.0),
        -Vector::Y,
        Vector::new(-0.3, -1.0, -0.1),
        Vector::new(-1.0, -0.5, -0.2),
        Vector::new(0.0, -0.5, -1.0),
        Vector::new(-0.5, 0.4, -1.0),
    ]
    .map(|dir| pos1.rotation * dir.normalize());
    let options = [
        options(0.31, 0.0, false),
        options(3.0, 0.0, true),
        options(3.0, 0.05, false),
        options(100.0, 0.0, true),
    ];
    let shapes = shapes();

    let (mut hits, mut misses) = (0, 0);
    for ground in &grounds {
        for (k, shape) in shapes.iter().enumerate() {
            for x in [-4, -1, 2] {
                for z in [-3, 2] {
                    for lift in [0.5, 0.95, 1.8] {
                        let origin = pos1 * Vector::new(x as Real + 0.37, lift, z as Real + 0.21);
                        let rotation = Rotation::from_rotation_x(0.1 * (x + k as i32) as Real);
                        let pos2 = Pose::from_parts(origin, rotation);
                        for dir in directions {
                            for options in options {
                                if assert_nearest_triangle_hit(
                                    &pos1,
                                    ground,
                                    &pos2,
                                    dir * 2.0,
                                    &**shape,
                                    options,
                                ) {
                                    hits += 1;
                                } else {
                                    misses += 1;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    assert!(
        hits > 1_000 && misses > 1_000,
        "{hits} hits, {misses} misses"
    );
}

/// A box some 100 across, cast at a target distance of about 9 over a cell 200 wide, at a
/// time of impact where GJK leaves the box about 8e-3 farther than the target distance from one
/// triangle: GJK's tolerance grows with the target distance and the length of the cast.
#[test]
fn heightfield_shape_cast_gives_the_nearest_triangle_hit_at_a_large_scale() {
    let heights = vec![-51.100224, 14.935101, -0.4987359, -7.1095943];
    let cell = HeightField::new(Array2::new(2, 2, heights), Vector::new(200.0, 1.0, 200.0));
    let cuboid = Cuboid::new(Vector::new(13.917611, 40.03146, 25.72583));
    let pos2 = Pose::from_parts(
        Vector::new(34.971046, 98.25103, -77.62439),
        Rotation::from_xyzw(-0.28899258, -0.5181429, -0.7841668, 0.18191668),
    );
    let vel2 = Vector::new(68.49045, -97.044716, -51.755314);

    let hits = assert_nearest_triangle_hit(
        &Pose::IDENTITY,
        &cell,
        &pos2,
        vel2,
        &cuboid,
        options(3.7495742, 9.270905, true),
    );

    assert!(hits);
}

/// Parry's dispatcher, recording each triangle it shape-casts against.
#[derive(Default)]
struct TriangleRecorder(Mutex<Vec<Triangle>>);

impl QueryDispatcher for TriangleRecorder {
    fn intersection_test(
        &self,
        pos12: &Pose,
        g1: &dyn Shape,
        g2: &dyn Shape,
    ) -> Result<ShapeIntersection, Unsupported> {
        DefaultQueryDispatcher.intersection_test(pos12, g1, g2)
    }

    fn distance(
        &self,
        pos12: &Pose,
        g1: &dyn Shape,
        g2: &dyn Shape,
    ) -> Result<ShapeDistance, Unsupported> {
        DefaultQueryDispatcher.distance(pos12, g1, g2)
    }

    fn contact(
        &self,
        pos12: &Pose,
        g1: &dyn Shape,
        g2: &dyn Shape,
        prediction: Real,
    ) -> Result<Option<Contact>, Unsupported> {
        DefaultQueryDispatcher.contact(pos12, g1, g2, prediction)
    }

    fn closest_points(
        &self,
        pos12: &Pose,
        g1: &dyn Shape,
        g2: &dyn Shape,
        max_dist: Real,
    ) -> Result<ClosestPoints, Unsupported> {
        DefaultQueryDispatcher.closest_points(pos12, g1, g2, max_dist)
    }

    fn cast_shapes(
        &self,
        pos12: &Pose,
        local_vel12: Vector,
        g1: &dyn Shape,
        g2: &dyn Shape,
        options: ShapeCastOptions,
    ) -> Result<Option<ShapeCastHit>, Unsupported> {
        if let Some(tri) = g1.as_triangle() {
            self.0.lock().unwrap().push(*tri);
        }
        DefaultQueryDispatcher.cast_shapes(pos12, local_vel12, g1, g2, options)
    }

    fn cast_shapes_nonlinear(
        &self,
        motion1: &NonlinearRigidMotion,
        g1: &dyn Shape,
        motion2: &NonlinearRigidMotion,
        g2: &dyn Shape,
        start_time: Real,
        end_time: Real,
        stop_at_penetration: bool,
    ) -> Result<Option<ShapeCastHit>, Unsupported> {
        DefaultQueryDispatcher.cast_shapes_nonlinear(
            motion1,
            g1,
            motion2,
            g2,
            start_time,
            end_time,
            stop_at_penetration,
        )
    }
}

/// The triangles of `heightfield` that a shape-cast of `g2` from `pos12` along `vel12` shape-casts
/// against, and its result.
fn cast_triangles(
    heightfield: &HeightField,
    pos12: &Pose,
    vel12: Vector,
    g2: &dyn Shape,
    options: ShapeCastOptions,
) -> (Vec<Triangle>, Option<ShapeCastHit>) {
    let recorder = TriangleRecorder::default();
    let hit =
        cast_shapes_heightfield_shape(&recorder, pos12, vel12, heightfield, g2, options).unwrap();
    (recorder.0.into_inner().unwrap(), hit)
}

/// A ball dropping by 0.1 toward the bottom of a valley, from 0.2 above it.
#[test]
fn heightfield_shape_cast_skips_the_triangles_out_of_reach() {
    let valley = heightfield(41, 40.0, |i, _| 0.02 * (i - 20.0).powi(2));
    let ball = Ball::new(0.4);
    let pos12 = Pose::translation(0.3, 0.2 + 0.4, 0.0);

    let (triangles, hit) = cast_triangles(
        &valley,
        &pos12,
        -Vector::Y,
        &ball,
        ShapeCastOptions::with_max_time_of_impact(0.1),
    );

    assert!(hit.is_none());
    assert!(triangles.is_empty(), "{} triangles cast", triangles.len());
}

/// A ball skimming down onto flat ground along -X hits it 10 along X from its start, at X = 5.
#[test]
fn heightfield_shape_cast_skips_the_triangles_past_the_nearest_hit() {
    let flat = heightfield(41, 40.0, |_, _| 0.0);
    let ball = Ball::new(0.4);
    let pos12 = Pose::translation(15.0, 1.4, 0.3);
    let vel12 = Vector::new(-1.0, -0.1, 0.0);

    let (triangles, hit) = cast_triangles(
        &flat,
        &pos12,
        vel12,
        &ball,
        ShapeCastOptions::with_max_time_of_impact(30.0),
    );

    let hit = hit.expect("the ball lands on the ground");
    assert_relative_eq!(hit.time_of_impact, 10.0, epsilon = 1.0e-3);
    // The ball's Aabb reaches X = 4.6 by then, short of the cell [3, 4].
    let past_the_hit = triangles
        .iter()
        .filter(|tri| tri.a.x.max(tri.b.x).max(tri.c.x) <= 4.0)
        .count();
    assert_eq!(past_the_hit, 0, "{} triangles cast", triangles.len());
}

/// A ball dropping onto flat ground lands at the very end of its shape-cast: the cast's
/// `max_time_of_impact` is the time of impact against the triangles on their own.
#[test]
fn heightfield_shape_cast_hits_at_the_end_of_the_cast() {
    let flat = heightfield(11, 10.0, |_, _| 0.0);
    let ball = Ball::new(0.4);
    let pos2 = Pose::translation(0.3, 1.4, 0.2);
    let landing = triangle_hits(
        &Pose::IDENTITY,
        &flat,
        &pos2,
        -Vector::Y,
        &ball,
        ShapeCastOptions::with_max_time_of_impact(2.0),
    )
    .iter()
    .map(|hit| hit.time_of_impact)
    .reduce(Real::min)
    .expect("the ball lands on the ground");

    let hits = assert_nearest_triangle_hit(
        &Pose::IDENTITY,
        &flat,
        &pos2,
        -Vector::Y,
        &ball,
        ShapeCastOptions::with_max_time_of_impact(landing),
    );

    assert!(hits);
}
