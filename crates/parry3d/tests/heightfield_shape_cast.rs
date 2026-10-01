// Shape-casting against a `HeightField` only shape-casts the triangles whose Aabb can be
// reached by the moving shape's Aabb. These check that no triangle that could be hit is
// skipped, and measure the gain.

use parry3d::math::{Pose, Real, Vector};
use parry3d::query::{cast_shapes, ShapeCastHit, ShapeCastOptions};
use parry3d::shape::{Ball, Capsule, Cuboid, HeightField, Shape};
use parry3d::utils::Array2;

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

    fn vector(&mut self, min: f32, max: f32) -> Vector {
        Vector::new(
            self.real(min, max),
            self.real(min, max),
            self.real(min, max),
        )
    }
}

/// Shape-casts every triangle of the heightfield, in the order the heightfield visits them.
fn brute_force_cast(
    heightfield: &HeightField,
    pos2: &Pose,
    vel2: Vector,
    g2: &dyn Shape,
    options: ShapeCastOptions,
) -> Option<ShapeCastHit> {
    let mut best_hit = None::<ShapeCastHit>;
    for tri in heightfield.triangles() {
        if let Some(hit) = cast_shapes(
            &Pose::identity(),
            Vector::ZERO,
            &tri,
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
fn heightfield_shape_cast_skips_no_reachable_triangle() {
    // Single-cell heightfields, so the heightfield's own traversal visits both triangles and
    // any difference with brute force comes from a skipped triangle.
    let mut rng = Lcg(0x5ca57);
    let shapes: [Box<dyn Shape>; 3] = [
        Box::new(Ball::new(0.4)),
        Box::new(Cuboid::new(Vector::new(0.6, 0.3, 0.9))),
        Box::new(Capsule::new_y(0.5, 0.3)),
    ];
    let mut num_hits = 0;

    for i in 0..20_000 {
        let heights = (0..4).map(|_| rng.real(-1.0, 1.0)).collect();
        let width = rng.real(0.5, 4.0);
        let heightfield = HeightField::new(
            Array2::new(2, 2, heights),
            Vector::new(width, rng.real(0.1, 2.0), width),
        );
        let g2 = &*shapes[i % shapes.len()];
        let pos2 = Pose::new(
            Vector::new(
                rng.real(-3.0, 3.0),
                rng.real(-2.0, 3.0),
                rng.real(-3.0, 3.0),
            ),
            rng.vector(-3.0, 3.0),
        );
        let vel2 = match i % 4 {
            // Mostly horizontal, like a body skimming over the terrain.
            0 => Vector::new(
                rng.real(-20.0, 20.0),
                rng.real(-1.0, 1.0),
                rng.real(-20.0, 20.0),
            ),
            // Straight down.
            1 => Vector::new(0.0, rng.real(-10.0, 0.0), 0.0),
            _ => rng.vector(-10.0, 10.0),
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

#[test]
#[ignore = "benchmark, run manually with --ignored --nocapture --test-threads=1"]
fn heightfield_shape_cast_perf() {
    // A car-sized box riding 20 cm above rolling terrain with 12.5 cm cells and sweeping one
    // 60 Hz step at 30 m/s (a CCD sweep that rarely hits). The terrain's Aabb contains the box,
    // so every triangle under it used to be shape-cast. Then balls falling onto the terrain,
    // where every triangle under the ball can be hit and none is skipped.
    let mut rng = Lcg(0xbe7c4);
    let n = 256;
    let size = 32.0;
    let cell = size / n as f32;
    let terrain_height = |x: f32, z: f32| 1.5 * (x / 5.0).sin() * (z / 7.0).cos();
    let mut heights = Array2::zeros(n + 1, n + 1);
    for i in 0..=n {
        for j in 0..=n {
            // Rows go along z, columns along x.
            let x = -size / 2.0 + j as f32 * cell;
            let z = -size / 2.0 + i as f32 * cell;
            heights[(i, j)] = terrain_height(x, z) + rng.real(-0.02, 0.02);
        }
    }
    let heightfield = HeightField::new(heights, Vector::new(size, 1.0, size));
    let chassis = Cuboid::new(Vector::new(0.9, 0.35, 2.2));
    let ball = Ball::new(0.3);

    let sweeps: Vec<(Pose, Vector)> = (0..10_000)
        .map(|_| {
            let mut pos = Pose::new(
                Vector::new(rng.real(-12.0, 12.0), 0.0, rng.real(-12.0, 12.0)),
                Vector::new(0.0, rng.real(-3.14, 3.14), 0.0),
            );
            // Lift the box 20 cm above the highest terrain under it.
            let aabb = chassis.compute_aabb(&pos);
            let mut ground: f32 = -100.0;
            let mut x = aabb.mins.x - cell;
            while x <= aabb.maxs.x + cell {
                let mut z = aabb.mins.z - cell;
                while z <= aabb.maxs.z + cell {
                    ground = ground.max(terrain_height(x, z) + 0.02);
                    z += cell / 2.0;
                }
                x += cell / 2.0;
            }
            pos.translation.y = ground + 0.2 + chassis.half_extents.y;
            let angle = rng.real(-3.14, 3.14);
            let vel = Vector::new(30.0 * angle.cos(), rng.real(-1.0, 1.0), 30.0 * angle.sin());
            (pos, vel)
        })
        .collect();
    let drops: Vec<Pose> = (0..10_000)
        .map(|_| {
            Pose::from_translation(Vector::new(
                rng.real(-12.0, 12.0),
                rng.real(2.0, 5.0),
                rng.real(-12.0, 12.0),
            ))
        })
        .collect();

    let time = |casts: &mut dyn FnMut() -> usize| {
        let t = std::time::Instant::now();
        let num_hits = casts();
        (t.elapsed().as_secs_f64() * 1.0e6 / 10_000.0, num_hits)
    };

    let (t_sweeps, sweep_hits) = time(&mut || {
        sweeps
            .iter()
            .filter(|(pos, vel)| {
                cast_shapes(
                    &Pose::identity(),
                    Vector::ZERO,
                    &heightfield,
                    pos,
                    *vel,
                    &chassis,
                    ShapeCastOptions::with_max_time_of_impact(1.0 / 60.0),
                )
                .unwrap()
                .is_some()
            })
            .count()
    });
    let (t_drops, drop_hits) = time(&mut || {
        drops
            .iter()
            .filter(|pos| {
                cast_shapes(
                    &Pose::identity(),
                    Vector::ZERO,
                    &heightfield,
                    pos,
                    Vector::new(0.0, -10.0, 0.0),
                    &ball,
                    ShapeCastOptions::with_max_time_of_impact(1.0),
                )
                .unwrap()
                .is_some()
            })
            .count()
    });

    println!(
        "chassis sweeps: {t_sweeps:.3}us/cast ({sweep_hits} hits), \
         falling balls: {t_drops:.3}us/cast ({drop_hits} hits)"
    );
}
