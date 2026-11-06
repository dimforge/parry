// Shape-casting against a `HeightField` skips the segments the moving shape can't hit within
// `max_time_of_impact`, or before the nearest hit found so far. These check that the result is
// still the nearest of the hits of a shape-cast against each segment on its own.

use approx::assert_relative_eq;
use parry2d::math::{Pose, Real, Vector};
use parry2d::query::details::cast_shapes_heightfield_shape;
use parry2d::query::{
    self, ClosestPoints, Contact, DefaultQueryDispatcher, NonlinearRigidMotion, QueryDispatcher,
    ShapeCastHit, ShapeCastOptions, ShapeDistance, ShapeIntersection, Unsupported,
};
use parry2d::shape::{Ball, Capsule, Cuboid, HeightField, RoundCuboid, Segment, Shape};
use std::sync::Mutex;

/// A heightfield `size` wide with `samples` samples, with the given height function.
fn heightfield(samples: usize, size: Real, height: impl Fn(Real) -> Real) -> HeightField {
    let heights = (0..samples).map(|i| height(i as Real)).collect();
    HeightField::new(heights, Vector::new(size, 1.0))
}

/// The shape-casts against each segment of `heightfield` on its own that hit.
fn segment_hits(
    pos1: &Pose,
    heightfield: &HeightField,
    pos2: &Pose,
    vel2: Vector,
    g2: &dyn Shape,
    options: ShapeCastOptions,
) -> Vec<ShapeCastHit> {
    heightfield
        .segments()
        .filter_map(|seg| {
            query::cast_shapes(pos1, Vector::ZERO, &seg, pos2, vel2, g2, options).unwrap()
        })
        .collect()
}

/// Shape-casts `g2` from `pos2` along `vel2` against `heightfield` at `pos1`, asserts that the
/// result is the nearest hit of the shape-casts against each segment on its own (any one among
/// exact ties), and returns whether it hits.
fn assert_nearest_segment_hit(
    pos1: &Pose,
    heightfield: &HeightField,
    pos2: &Pose,
    vel2: Vector,
    g2: &dyn Shape,
    options: ShapeCastOptions,
) -> bool {
    let hit = query::cast_shapes(pos1, Vector::ZERO, heightfield, pos2, vel2, g2, options).unwrap();
    let hits = segment_hits(pos1, heightfield, pos2, vel2, g2, options);
    let nearest = hits.iter().map(|hit| hit.time_of_impact).reduce(Real::min);
    let context = || format!("{pos2:?}, {vel2:?}, {options:?}");

    assert_eq!(hit.map(|hit| hit.time_of_impact), nearest, "{}", context());
    if let Some(hit) = hit {
        // `ShapeCastHit` has no `PartialEq`; `Debug` prints each float so that it reads back
        // exactly.
        let hit = format!("{hit:?}");
        assert!(
            hits.iter().any(|seg_hit| format!("{seg_hit:?}") == hit),
            "{hit} is no segment's hit: {}",
            context()
        );
    }
    hit.is_some()
}

fn shapes() -> [Box<dyn Shape>; 4] {
    [
        Box::new(Ball::new(0.4)),
        Box::new(Cuboid::new(Vector::new(0.6, 0.3))),
        Box::new(Capsule::new_y(0.5, 0.4)),
        Box::new(RoundCuboid {
            inner_shape: Cuboid::new(Vector::new(0.3, 0.2)),
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
fn heightfield_shape_cast_gives_the_nearest_segment_hit() {
    // Waves about 10 m long and 0.5 m high, and ground that rises and falls by up to about 1 m
    // from one sample to the next.
    let grounds = [
        heightfield(41, 40.0, |i| 0.5 * (0.6 * i).cos()),
        heightfield(41, 40.0, |i| 0.6 * ((1.3 * i).sin() + (0.7 * i).cos())),
    ];
    // The heightfield is posed off the origin and turned, so that the cast runs in its frame.
    let pos1 = Pose::new(Vector::new(3.0, -1.0), 0.3);
    // In the heightfield's frame.
    let directions = [
        Vector::new(1.0, 0.0),
        Vector::new(-1.0, -0.01),
        Vector::new(1.0, 0.02),
        -Vector::Y,
        Vector::new(0.3, -1.0),
        Vector::new(-1.0, -1.0),
        Vector::new(1.0, -0.5),
        Vector::new(-0.5, 0.4),
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
            for x in (-19..=19).step_by(2) {
                for lift in [0.5, 0.95, 1.8] {
                    let origin = pos1 * Vector::new(x as Real + 0.37, lift);
                    let pos2 = Pose::new(origin, 0.1 * (x + k as i32) as Real);
                    for dir in directions {
                        for options in options {
                            if assert_nearest_segment_hit(
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
    assert!(
        hits > 1_000 && misses > 1_000,
        "{hits} hits, {misses} misses"
    );
}

/// Parry's dispatcher, recording each segment it shape-casts against.
#[derive(Default)]
struct SegmentRecorder(Mutex<Vec<Segment>>);

impl QueryDispatcher for SegmentRecorder {
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
        if let Some(seg) = g1.as_segment() {
            self.0.lock().unwrap().push(*seg);
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

/// The segments of `heightfield` that a shape-cast of `g2` from `pos12` along `vel12` shape-casts
/// against, and its result.
fn cast_segments(
    heightfield: &HeightField,
    pos12: &Pose,
    vel12: Vector,
    g2: &dyn Shape,
    options: ShapeCastOptions,
) -> (Vec<Segment>, Option<ShapeCastHit>) {
    let recorder = SegmentRecorder::default();
    let hit =
        cast_shapes_heightfield_shape(&recorder, pos12, vel12, heightfield, g2, options).unwrap();
    (recorder.0.into_inner().unwrap(), hit)
}

/// A ball dropping by 0.1 toward the bottom of a valley, from 0.2 above it.
#[test]
fn heightfield_shape_cast_skips_the_segments_out_of_reach() {
    let valley = heightfield(41, 40.0, |i| 0.02 * (i - 20.0).powi(2));
    let ball = Ball::new(0.4);
    let pos12 = Pose::translation(0.3, 0.2 + 0.4);

    let (segments, hit) = cast_segments(
        &valley,
        &pos12,
        -Vector::Y,
        &ball,
        ShapeCastOptions::with_max_time_of_impact(0.1),
    );

    assert!(hit.is_none());
    assert!(segments.is_empty(), "{} segments cast", segments.len());
}

/// A ball skimming down onto flat ground along -X hits it 10 along X from its start, at X = 5.
#[test]
fn heightfield_shape_cast_skips_the_segments_past_the_nearest_hit() {
    let flat = heightfield(41, 40.0, |_| 0.0);
    let ball = Ball::new(0.4);
    let pos12 = Pose::translation(15.0, 1.4);
    let vel12 = Vector::new(-1.0, -0.1);

    let (segments, hit) = cast_segments(
        &flat,
        &pos12,
        vel12,
        &ball,
        ShapeCastOptions::with_max_time_of_impact(30.0),
    );

    let hit = hit.expect("the ball lands on the ground");
    assert_relative_eq!(hit.time_of_impact, 10.0, epsilon = 1.0e-3);
    // The ball's Aabb reaches X = 4.6 by then, short of the segment [3, 4].
    let past_the_hit = segments
        .iter()
        .filter(|seg| seg.a.x.max(seg.b.x) <= 4.0)
        .count();
    assert_eq!(past_the_hit, 0, "{} segments cast", segments.len());
}

/// A ball dropping onto flat ground lands at the very end of its shape-cast: the cast's
/// `max_time_of_impact` is the time of impact against the segments on their own.
#[test]
fn heightfield_shape_cast_hits_at_the_end_of_the_cast() {
    let flat = heightfield(11, 10.0, |_| 0.0);
    let ball = Ball::new(0.4);
    let pos2 = Pose::translation(0.3, 1.4);
    let landing = segment_hits(
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

    let hits = assert_nearest_segment_hit(
        &Pose::IDENTITY,
        &flat,
        &pos2,
        -Vector::Y,
        &ball,
        ShapeCastOptions::with_max_time_of_impact(landing),
    );

    assert!(hits);
}

/// Shapes standing on flat ground, or a hair above it, where GJK may report a hit on a segment the
/// shape stays clear of, slide along it.
#[test]
fn heightfield_shape_cast_gives_the_nearest_segment_hit_within_gjks_tolerance() {
    let flat = heightfield(11, 10.0, |_| 0.0);
    let shapes = shapes();
    let half_heights = [0.4, 0.3, 0.9, 0.35];
    let directions = [
        Vector::X,
        -Vector::X,
        Vector::new(1.0, 0.01),
        Vector::new(-1.0, -0.001),
    ];

    for (shape, half_height) in shapes.iter().zip(half_heights) {
        for hair in [0.0, 1.0e-6, 1.0e-4] {
            let pos2 = Pose::translation(0.3, half_height + hair);
            for dir in directions {
                for options in [options(1.0, 0.0, true), options(1.0, 0.0, false)] {
                    assert_nearest_segment_hit(
                        &Pose::IDENTITY,
                        &flat,
                        &pos2,
                        dir.normalize(),
                        &**shape,
                        options,
                    );
                }
            }
        }
    }
}
