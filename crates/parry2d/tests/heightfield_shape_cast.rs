// Shape-casting against a `HeightField` only shape-casts the segments whose Aabb can be
// reached by the moving shape's Aabb. This checks that no segment that could be hit is
// skipped.

use parry2d::math::{Pose, Real, Vector};
use parry2d::query::{cast_shapes, ShapeCastHit, ShapeCastOptions};
use parry2d::shape::{Ball, Capsule, Cuboid, HeightField, Shape};

/// Tiny deterministic LCG so the tests don't depend on rand seeding.
struct Lcg(u64);

impl Lcg {
    fn next_u32(&mut self) -> u32 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (self.0 >> 32) as u32
    }

    fn real(&mut self, min: f32, max: f32) -> f32 {
        min + (max - min) * (self.next_u32() as f32 / u32::MAX as f32)
    }
}

/// Shape-casts every segment of the heightfield, in the order of the heightfield.
fn brute_force_cast(
    heightfield: &HeightField,
    pos2: &Pose,
    vel2: Vector,
    g2: &dyn Shape,
    options: ShapeCastOptions,
) -> Option<ShapeCastHit> {
    let mut best_hit = None::<ShapeCastHit>;
    for seg in heightfield.segments() {
        if let Some(hit) = cast_shapes(
            &Pose::identity(),
            Vector::ZERO,
            &seg,
            pos2,
            vel2,
            g2,
            options,
        )
        .unwrap()
        {
            if hit.time_of_impact < best_hit.map_or(Real::MAX, |h| h.time_of_impact) {
                best_hit = Some(hit);
            }
        }
    }
    best_hit
}

#[test]
fn heightfield_shape_cast_skips_no_reachable_segment() {
    // Single-segment heightfields with the shape starting over the segment, so the
    // heightfield's own traversal visits the segment and any difference with brute force
    // comes from a skipped segment.
    let mut rng = Lcg(0x5ca57);
    let shapes: [Box<dyn Shape>; 3] = [
        Box::new(Ball::new(0.4)),
        Box::new(Cuboid::new(Vector::new(0.6, 0.3))),
        Box::new(Capsule::new_y(0.5, 0.3)),
    ];
    let mut num_hits = 0;

    for i in 0..20_000 {
        let heights = (0..2).map(|_| rng.real(-1.0, 1.0)).collect();
        let width = rng.real(0.5, 4.0);
        let heightfield = HeightField::new(heights, Vector::new(width, rng.real(0.1, 2.0)));
        let g2 = &*shapes[i % shapes.len()];
        let pos2 = Pose::new(
            Vector::new(rng.real(-width, width) / 2.0, rng.real(-2.0, 3.0)),
            rng.real(-3.0, 3.0),
        );
        let vel2 = match i % 4 {
            // Mostly horizontal, like a body skimming over the terrain.
            0 => Vector::new(rng.real(-20.0, 20.0), rng.real(-1.0, 1.0)),
            // Straight down.
            1 => Vector::new(0.0, rng.real(-10.0, 0.0)),
            _ => Vector::new(rng.real(-10.0, 10.0), rng.real(-10.0, 10.0)),
        };
        let options = ShapeCastOptions {
            max_time_of_impact: if i % 5 == 0 {
                Real::MAX
            } else {
                rng.real(0.01, 2.0)
            },
            target_distance: if i % 3 == 0 { rng.real(0.0, 0.2) } else { 0.0 },
            stop_at_penetration: i % 7 != 0,
            compute_impact_geometry_on_penetration: true,
        };

        let hit = cast_shapes(
            &Pose::identity(),
            Vector::ZERO,
            &heightfield,
            &pos2,
            vel2,
            g2,
            options,
        )
        .unwrap();
        let brute_hit = brute_force_cast(&heightfield, &pos2, vel2, g2, options);

        assert_eq!(
            hit.map(|h| (h.time_of_impact, h.witness1, h.normal1)),
            brute_hit.map(|h| (h.time_of_impact, h.witness1, h.normal1)),
            "cast {i}: {pos2:?}, {vel2:?}, {options:?}"
        );
        num_hits += hit.is_some() as usize;
    }

    // Make sure the casts above are neither all misses nor all hits.
    assert!(num_hits > 2_000 && num_hits < 18_000, "{num_hits} hits");
}
